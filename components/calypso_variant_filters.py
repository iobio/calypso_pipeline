#!/usr/bin/python

from __future__ import print_function
from datetime import date
import json
import math
import os

# Process all the variant filters in Calypso and add to the project
def variantFilters(mosaicConfig, rootPath, samples, privateAnnotations, projectId, api_vf):
  existingFilters = {}

  # Define the path to the variant filter json file
  jsonPath = rootPath + str('/mosaic_variant_filters/')

  # Get all of the filters that exist in the project
  try:
    for existingFilter in json.loads(os.popen(api_vf.getVariantFilters(mosaicConfig, projectId)).read()): existingFilters[existingFilter['name']] = existingFilter['id']
  except: fail('Unable to get existing variant filters for project ' + str(projectId))

  # Get the ids of the samples
  probandId, motherId, fatherId = getSampleInfo(samples)

  # Loop over all the filter json files
  for filterJson in os.listdir(jsonPath):
    if filterJson.endswith('.json'):
      filterData = openFilterFile(jsonPath, filterJson)
      filterName = getFilterName(filterData, filterJson)

      # The 'filters' field needs to be present in the json
      if 'filters' not in filterData: fail('Variant filter json file with the name ' + str(filterJson) + ' does not contain the required \'filters\' section')
      if filterName in existingFilters: deleteFilter(mosaicConfig, filterName, projectId, existingFilters[filterName], api_vf)
      filterData = getGenotypeFilter(filterData, filterJson, probandId, motherId, fatherId)
      filterData = annotationFilters(filterData, filterName, privateAnnotations)

      # Create the annotation filter
      createFilter(filterData, mosaicConfig, projectId, api_vf)

# Get information about the samples
def getSampleInfo(samples):
  probandId = False
  motherId  = False
  fatherId  = False

  # Get the id of the proband and the parents
  for sample in samples:
    if samples[sample]["relationship"] == "Proband": probandId = samples[sample]["mosaicId"]
    elif samples[sample]["relationship"] == "Mother": motherId = samples[sample]["mosaicId"]
    elif samples[sample]["relationship"] == "Father": fatherId = samples[sample]["mosaicId"]

  # Return the ids
  return probandId, motherId, fatherId

# Open the annotation filter file
def openFilterFile(filepath, filename):

  # Get the information from the json file defining the filter
  try: filterFile = open(filepath + filename, 'r')
  except: fail('File ' + filepath + filename + ' could not be opened')
  try: data = json.load(filterFile)
  except: fail('File ' + filepath + filename + ' is not a valid json file')
  filterFile.close()

  # Return the annotation filter information
  return data

# Get the name of the filter to add
def getFilterName(data, filename):

  # Get the name of the filter to create
  try: filterName = data['name']
  except: fail('Variant filter file, ' + str(filename) + ', defines a filter with no name. A name needs to be supplied')

  # Return the name of the filter
  return filterName

# Delete an existing filter
def deleteFilter(mosaicConfig, filterName, projectId, filterId, api_vf):
  try: deleteFilter = os.popen(api_vf.deleteVariantFilter(mosaicConfig, projectId, filterId)).read()
  except: fail('Unable to delete variant filter, ' + str(filterName) + ' from project ' + str(projectId))

# Get information on the genotype filters
def getGenotypeFilter(data, filename, probandId, motherId, fatherId):

  # Store the allowed genotype options for saved filters
  genotypeOptions = []
  genotypeOptions.append('ref_samples')
  genotypeOptions.append('alt_samples')
  genotypeOptions.append('het_samples')
  genotypeOptions.append('hom_samples')

  # Check what genotype filters need to be applied
  try: genotypes = data['genotypes']
  except: fail('Mosaic variant filter json file, ' + str(filterJson) + ', does not contain the "genotypes" field')
  for genotype in genotypes:
    if genotype not in genotypeOptions: fail('Mosaic variant filter json file, ' + str(filename) + ', contains unknown genotype options: ' + str(genotype))
    if not genotypes[genotype]: continue

    # Check which samples need to have the requested genotype and add to the command
    sampleList = []
    for sample in genotypes[genotype]:
      if sample == 'proband': sampleList.append(probandId)
      elif sample == 'mother':
        if motherId: sampleList.append(motherId)
      elif sample == 'father':
        if fatherId: sampleList.append(fatherId)
      elif sample == 'sibling':
        fail('Can\'t handle siblings in filters yet')
      else: fail('Unknown sample, ' + str(sample) + ' in genotypes for Mosaic variant filter (' + str(filename) + ')')

    # Add the genotype filter to the filters listed in the json
    data['filters'][genotype] = sampleList

  # Return the updated data
  return data

# Process the annotation filters
def annotationFilters(data, filename, privateAnnotations):
  removeAnnotations = []

  # Make sure the annotation_filters section exists
  if 'annotation_filters' not in data['filters']: fail('Annotation filter json file with the name ' + str(filename) + ' does not contain the required \'annotation_filters\' section')

  # Check the filters provided in the json. The annotation filters need to extract the
  # uids for the annotations
  for index, annotationFilter in enumerate(data['filters']['annotation_filters']):

    # If the uid is provided with the filter, this annotation is complete. If not, the name
    # field must be present, and this must point to a created annotation uid. Remove the name
    # field and replace it with the uid.
    if 'uid' not in annotationFilter:
      try: annotationName = annotationFilter['name']
      except: fail('Annotation ' + str(filterJson) + ' does not have a uid, but also no name')
      annotationFilter.pop('name')

      # If the annotation is listed as optional and the name doesn't point to a created annotation
      # delete this annotation from the filter. This could be, for example, a filter on genotype
      # quality for the mother, but the mother isn't present in the project
      isOptional = annotationFilter.pop('optional') if 'optional' in annotationFilter else False

      # Check if the annotation has been created in this execution of Calypso
      foundMatch = False
      for annotation in privateAnnotations:
        if privateAnnotations[annotation]['name'] == annotationName:
          annotationFilter['uid'] = annotation
          foundMatch = True
          break

      # If no annotation with the required name was created and the annotation is optional, remove the annotation
      if not foundMatch:
        if isOptional: removeAnnotations.append(index)
        else: fail('Annotation, ' + str(annotationName) + ', in ' + str(filename) + ' does not have a uid, and the name does not point to a created annotation')

  # Remove any optional filters that did not have uids
  for index in reversed(removeAnnotations): data['filters']['annotation_filters'].pop(index)

  # Return the updated annotation filter data
  return data

# Create the annotation filter
def createFilter(data, mosaicConfig, projectId, api_vf):
  try: execute = json.loads(os.popen(api_vf.postVariantFilter(mosaicConfig, data['name'], data['filters'], projectId)).read())
  except: fail('Failed to create annotation filter ' + data['name'])

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
