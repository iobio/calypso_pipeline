from __future__ import print_function
from os.path import exists
from datetime import date
import json
import math
import os

# Process the json file describing the filters to apply
def readRequiredFilters(mosaicConfig, mosaicInfo, dataDirectory, samples, projectAnnotations, projectId, api_ps, api_vf):
  categories = {}
  filters    = {}

  # Create a mapping from the samples relationship (e.g. proband, or mother), to a Mosaic sample id
  sampleMap = {}
  for sample in samples: sampleMap[samples[sample]['relationship'].lower()] = samples[sample]['mosaicId']

  # We also need an annotation map that links the current version of ClinVar to a uid for the filter
  annotationMap = {}
  annotationMap['clinvar_latest'] = mosaicInfo['resources']['clinVar']['annotations']['CLNSIG']['uid']

  # Read the json
  filterInfo = readJson(str(dataDirectory) + str(mosaicInfo['variantFilters']))

  # The path to the files must be provided
  if 'path_to_files' not in filterInfo: fail('The json file describing variant filters is missing the "path_to_files" section')
  filterPath = str(filterInfo['path_to_files']) if filterInfo['path_to_files'].endswith('/') else str(filterInfo['path_to_files']) + '/'

  # Loop over the categories and validate all information
  if 'categories' not in filterInfo: fail('The json file describing variant filters is missing the "categories" section')
  for category in filterInfo['categories']:

    # For each category, loop over the names of the filters, and store their sort order. Check that there are no
    # duplicated sort positions
    categories[category] = {}
    for name in filterInfo['categories'][category]:
      position = filterInfo['categories'][category][name]
      if position in categories[category]: fail('Filter ' + str(name) + ' in category ' + str(category) + ' has the same sort position as a different filter')
      categories[category][position] = name
      if name in filters: fail('Filter "' + str(name) + '" appears multiple times in the filter description json')
      filters[name] = {'category': category, 'sortPosition': position}

  # Now check that all of the filters defined for each category are described in detail in the "filters" section of the json
  if 'filters' not in filterInfo: fail('The json file describing variant filters is missing the "filters" section')
  for category in categories:
    for position in categories[category]:
      name = categories[category][position]
      if name not in filterInfo['filters']: fail('Filter "' + str(name) + '" appears in filter category "' + str(category) + '", but is not described in the "filters" section')
    
      # Check the genotype information for the filter
      filters[name]['info'] = filterInfo['filters'][name]
      if 'genotypes' in filters[name]['info']: filters[name]['info'] = checkGenotypeFilters(filters[name]['info'], name, list(samples.keys()), sampleMap)
      filters[name]['info'] = checkAnnotationFilters(filters[name]['info'], name, projectAnnotations, annotationMap)

  # Get all of the filters that exist in the project, and check which of these share a name with a filter to be created
  deleteFilters(mosaicConfig, projectId, filters, api_vf)

  # Create all the required filters and update their categories and sort order in the project settings
  createFilters(mosaicConfig, projectId, categories, filters, api_ps, api_vf)

# Read the input json file to get all of the filters and their categories and sort orders
def readJson(filename):

  # Check that the file defining the filters exists
  if not exists(filename): fail('Could not find the json file ' + str(filename))

  # The file describing the variant filters should be in json format. Fail if the file is not valid
  try: jsonFile = open(filename, "r")
  except: fail('Could not open the json file: ' + str(filename))
  try: data = json.load(jsonFile)
  except: fail('Could not read contents of json file ' + str(filename) + '. Check that this is a valid json')

  # Return the json information
  return data

# Get information on the genotype filters
def checkGenotypeFilters(data, name, sampleIds, sampleMap):

  # Store the allowed genotype options for saved filters
  genotypeOptions = []
  genotypeOptions.append('ref_samples')
  genotypeOptions.append('alt_samples')
  genotypeOptions.append('het_samples')
  genotypeOptions.append('hom_samples')

  # Check what genotype filters need to be applied and that the supplied genotypes are valid
  for genotype in data['genotypes']:
    if genotype not in genotypeOptions: fail('Mosaic variant filter with the name ' + str(name) + ', contains an unknown genotype option: ' + str(genotype))
    if not data['genotypes'][genotype]: continue

    # Check which samples need to have the requested genotype and add to the command. Use the supplied sampleIds
    # list to check that these samples are in the project
    sampleList = []
    if type(data['genotypes'][genotype]) != list: fail('Mosaic variant filter with the name ' + str(name) + ' has an invalid genotypes section')
    for sample in data['genotypes'][genotype]:

      # The genotype filter must either contain a valid sample id for the project, or the value in the json (e.g. proband)
      # must be present in the sampleMap and point to a valid sample id for this project
      sampleId = sampleMap[sample] if sample in sampleMap else False
      if not sampleId:
        try: sampleId = int(sample)
        except: fail('Mosaic variant filter ' + str(name) + ' references a sample with a non-integer id: ' + str(sample))
        if int(sampleId) not in sampleIds: fail('Mosaic variant filter ' + str(name) + ' references sample ' + str(sample) + ' which is not in the requested project')
      sampleList.append(sampleId)

    # Add the genotype filter to the filters listed in the json
    data['filters'][genotype] = sampleList

  # Return the updated data
  return data

# Process the annotation filters
def checkAnnotationFilters(data, name, uids, annotationMap):
  annotationFilters = []

  # Make sure the annotation_filters section exists
  if 'annotation_filters' not in data['filters']: fail('Annotation filter ' + str(name) + ' does not contain the required "annotation_filters" section')

  # Check the filters provided in the json. The annotation filters need to extract the uids for the annotations, so
  # ensure that each annotation has a valid uid (e.g. it is present in the project), and that supporting information
  # e.g. a minimum value cannot be supplied for a string annotation, is valid
  for aFilter in data['filters']['annotation_filters']:

    # The json file must contain either a valid uid for a project annotation, the name of a valid private annotation (for
    # projects where private annotations are created, a new filter template shouldn't be required for every project), or
    # have a name in the annotation map to relate a name to a uid. This is used for annotations (e.g. ClinVar) that are
    # regularly updated, so the template does not need to be updated for updating annotations.
    # If a uid is provided, check it is valid
    uid = False
    if 'uid' in aFilter: uid = aFilter['uid']

    # If a name is provided instead of a uid...
    elif 'name' in aFilter:

      # ...check if this name is in the annotationMap and if so, use the mapped uid
      if aFilter['name'] in annotationMap:
        uid            = annotationMap[aFilter['name']]
        aFilter['uid'] = uid
        del aFilter['name']

      # ...or if the name is not in the annotationMap, check if a private annotation with this name exists in the project
      else:
        for pUid in uids:
          if str(uids[pUid]['name']) == str(aFilter['name']):
            uid            = pUid
            aFilter['uid'] = uid
            del aFilter['name']
            break

    # Check if this annotation is optional. If no uid has been found and the annotation is optional, it should be removed. This could
    # occur, for example, if the annotation is the genotype quality of the mother, but there is no mother for this case, and so the
    # private annotation was never created, and so no filter should be created for this
    isOptional = False
    if 'optional' in aFilter: isOptional = aFilter.pop('optional')
    if not uid and isOptional: continue
    elif not uid: fail('No uid can be determined for annotation filter ' + str(name))
    if 'include_nulls' not in aFilter: fail('Annotation filter ' + str(name) + ' contains a filter with no "include_nulls" section')

    # If the annotation is a string, the "values" field must be present
    if uids[uid]['type'] == 'string':
      if 'values' not in aFilter: fail('Annotation filter ' + str(name) + ' contains a string based filter with no "values" section')
      if type(aFilter['values']) != list: fail('Annotation filter ' + str(name) + ' contains a string based filter with a "values" section that is not a list')

    # If the annotation is a float, check that the specified operation is valid
    elif uids[uid] == 'float':

      # Loop over all the fields for the filter and check that they are valid
      hasRequiredValue = False
      for value in aFilter:
        if value == 'uid': continue
        elif value == 'include_nulls': continue

        # The filter can define a minimum value
        elif value == 'min':
          try: float(aFilter[value])
          except: fail('Annotation filter ' + str(name) + ' has a "min" that is not a float')
          hasRequiredValue = True

        # The filter can define a minimum value
        elif value == 'max':
          try: float(aFilter[value])
          except: fail('Annotation filter ' + str(name) + ' has a "max" that is not a float')
          hasRequiredValue = True

        # Other fields are not recognised
        else: fail('Annotation filter ' + str(name) + ' contains an unrecognised field: ' + str(value))

      # If no comparison fields were provided, fail
      if not hasRequiredValue: fail('Annotation filter ' + str(name) + ' contains a filter based on a float, but no comparison operators have been included')

    # Store the updated annotation filters
    annotationFilters.append(aFilter)

  # Replace the annotation filters with the updated version
  data['filters']['annotation_filters'] = annotationFilters

  # Return the updated annotation information
  return data

# Get all of the filters that exist in the project, and check which of these share a name with a filter to be created
def deleteFilters(mosaicConfig, projectId, filters, api_vf):
  existingFilters = api_vf.getVariantFilterIdsNames(mosaicConfig, projectId)
  removeFilters   = []
  for filterId in existingFilters: 
    if existingFilters[filterId] in filters.keys(): removeFilters.append(filterId)

  # Delete existing filters sharing a name with filters to be created
  for filterId in removeFilters: api_vf.deleteVariantFilter(mosaicConfig, projectId, filterId)

# Create all the required filters and update their categories and sort order in the project settings
def createFilters(mosaicConfig, projectId, categories, filters, api_ps, api_vf):
  sortedFilters = []
  for category in categories:
    record = {'category': category, 'sortOrder': []}
    for sortId in sorted(categories[category]):
      name = categories[category][sortId]

      # Create the filter
      filterId = api_vf.createVariantFilter(mosaicConfig, projectId, name, category, filters[name]['info']['filters'])
      record['sortOrder'].append(str(filterId))

    # Populate the object used to update the Mosaic project settings
    sortedFilters.append(record)

  # Set the sort orders for all the categories
  api_ps.setVariantFilterSortOrder(mosaicConfig, projectId, sortedFilters)

# Create the annotation filter
def createFilter(data, mosaicConfig, projectId, api_vf):
  try: execute = json.loads(os.popen(api_vf.postVariantFilter(mosaicConfig, data['name'], data['filters'], projectId)).read())
  except: fail('Failed to create annotation filter ' + data['name'])

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
