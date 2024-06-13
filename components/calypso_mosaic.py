from datetime import date
import json
import math
import os

# Get all of the project attributes in the target project
def getProjectAttributes(mosaicConfig, projectId, api_pa):
  projectAttributes = {}

  # Get the first page of attributes
  data = api_pa.getProjectAttributes(mosaicConfig, projectId)
  for attribute in data:

    # Get the value for this project
    for attributeProject in attribute['values']:
      value = False
      if int(attributeProject['project_id']) == int(projectId): value = attributeProject['value']
    projectAttributes[attribute['uid']] = {'name': attribute['name'], 'id': attribute['id'], 'value': value}

  # Return the project attributes
  return projectAttributes

# Get all of the public attributes that are available for import
def getPublicProjectAttributes(mosaicConfig, api_pa):
  publicAttributes = {}
  data = api_pa.getPublicProjectAttributesNameIdUid(mosaicConfig)
  for record in data: publicAttributes[record['uid']] = {'name': record['name'], 'id': record['id']}

  # Return the public project attributes
  return publicAttributes

# Get the resource version from the previous execution of the Calypso pipeline (if there is one)
def getPreviousResourceVersion(projectAttributes):
  version = False
  for attribute in projectAttributes:
    if projectAttributes[attribute]['name'] == 'Calypso resource version': version = projectAttributes[attribute]['value']

  # Return the version
  return version

# Get all available public annotations
def getPublicAnnotations(project, mosaicConfig, projectId, api_va):
  publicAnnotations = {}
  for annotation in project.get_variant_annotations_to_import():
    publicAnnotations[annotation['uid']] = {'name': annotation['name'], 'id': annotation['id'], 'type': annotation['value_type']}

  # Return the dictionary of available public annotations
  return publicAnnotations

# Create all the required private annotations
def createPrivateAnnotations(mosaicConfig, resources, projectAnnotations, samples, projectId, api_va):
  privateAnnotations = {}

  # Get the names (not the uids) of all the annotations in the project
  existingAnnotations = []
  for annotation in projectAnnotations: existingAnnotations.append(projectAnnotations[annotation]['name'])

  # Loop over each resource in the Mosaic resources file
  for resource in resources:
    annotations = {}

    # Loop over all annotations for this resource
    if resources[resource]['annotation_type'] == 'private':

      # The project id for these annotations is the target project as these annotations will be created there. Update
      # mosaicInfo to reflect this.
      resources[resource]['project_id'] = projectId

      # Loop over all annotations for this resource
      for annotation in resources[resource]['annotations']:
        valueType = resources[resource]['annotations'][annotation]['type']

        # If the annotation has "isSample" set, there needs to be an annotation for every sample in the vcf file for the 
        # project (some projects include samples that were not sequenced and so these will not have annotations). Append the
        # the samples relationship to the proband to the annotation name. 
        if resources[resource]['annotations'][annotation]['isSample']:
          for sample in samples:

            # Include the sample name in the annotation, otherwise non-unique names could result. For example, if the proband
            # has 2 brothers, the generated annotation name for the 2 brothers will be identical
            if samples[sample]['vcf_sample_name']:
              annotations[str(annotation) + ' ' + str(samples[sample]['relation']) + ' ' + str(sample)] = valueType
        else: annotations[annotation] = valueType

      # Remove this annotation from the resources dictionary. It will be replaced below with the new annotations with the
      # corresponding uids
      resources[resource]['annotations'] = {}

    # Loop over the annotations to be created, check that there doesn't already exist an annotation of that name in the
    # project, and if not, create a new, private annotation
    for annotation in annotations:
      valueType = annotations[annotation]
  
      # If the annotation already exists, get the uid of the annotation, otherwise create the annotation.
      if annotation in existingAnnotations:
        for existingAnnotation in projectAnnotations: 
          if str(projectAnnotations[existingAnnotation]['name']) == annotation:
            privateAnnotations[existingAnnotation] = {'name': annotation, 'id': projectAnnotations[existingAnnotation]['id']}
            resources[resource]['annotations'][annotation] = {'uid': existingAnnotation, 'type': valueType, 'id': projectAnnotations[existingAnnotation]['id']}
            break
      else:
        data = api_va.createPrivateAnnotationIdUid(mosaicConfig, annotation, valueType, projectId)
        privateAnnotations[data['uid']] = {'name': annotation, 'id': data['id']}
        resources[resource]['annotations'][annotation] = {'uid': data['uid'], 'type': valueType, 'id': data['id']}

        # Add the created private annotation to the projectAnnotations dictionary
        projectAnnotations[data['uid']] = {'id': data['id'], 'name': annotation, 'type': valueType}

  # Return the created private annotations
  return privateAnnotations, projectAnnotations

# Get all the project annotations already in a Mosaic project
def getProjectAnnotations(project, mosaicConfig, projectId, api_va):
  projectAnnotations = {}

  # Get the annotations
  data = api_va.getVariantAnnotationsNameIdUIdType(mosaicConfig, projectId)
  for annotation in data: projectAnnotations[annotation['uid']] = {'id': annotation['id'], 'name': annotation['name'], 'type': annotation['type']}

  # Return a formatted dictionary of annotations
  return projectAnnotations

# Remove any annotations from a Mosaic project. This must be provided with a list of uids
def removeAnnotations(mosaicConfig, uids, projectAnnotations, projectId, api_va):

  # Loop over the annotations to remove, get the annotation id, and remove them from the project
  for uid in uids:
    if uid in projectAnnotations: data = api_va.deleteAnnotationById(mosaicConfig, projectId, projectAnnotations[uid]['id'])

# Import a list of public annotations into a Mosaic project
def importAnnotations(mosaicConfig, resources, projectAnnotations, publicAnnotations, projectId, api_va):

  # Loop over all of the resources in the Mosaic resources file and get all of the annotation uids to import.
  # Only public attributes can be imported
  for resource in resources:
    if resources[resource]['annotation_type'] == 'public':

      # Loop over all the annotations required for this resource
      for annotation in resources[resource]['annotations']:
        uid = resources[resource]['annotations'][annotation]['uid']

        # Check if the annotation is already in the target project. If not, import the annotation and update the
        # projectAnnotations dictionary
        if uid not in projectAnnotations:
          api_va.importAnnotation(mosaicConfig, publicAnnotations[uid]['id'], projectId)
          projectAnnotations[uid] = {'id': publicAnnotations[uid]['id'], 'name': publicAnnotations[uid]['name'], 'type': publicAnnotations[uid]['type']}

  # Return the updated projectAnnotations dictionary
  return projectAnnotations

# Set the projects default annotations. The default annotations are what a new user will see in the table by default
# and all users will have the option to reset the annotations table to show only the default annotations
def defaultAnnotations(mosaicConfig, defaultAnnotations, projectAnnotations, publicAnnotations, privateAnnotations, projectId, api_ps):

  # Store all default annotation ids in a list
  annotationIds = []

  # Get the ids of all the default annotations
  for annotation in defaultAnnotations:
    if annotation in projectAnnotations:
      annotationIds.append(projectAnnotations[annotation]['id'])
    elif annotation in publicAnnotations:
      annotationIds.append(publicAnnotations[annotation]['id'])

    # Check if the annotation is a private annotation. In this case, the supplied annotation will not be the uid as
    # this is not known, but will be the annotation name
    else:
      pAnnoId = False
      for pAnno in privateAnnotations:
        if str(privateAnnotations[pAnno]['name']) == str(annotation): pAnnoId = privateAnnotations[pAnno]['id']
      if pAnnoId: annotationIds.append(pAnnoId)
      else: fail('Default annotation ' + str(annotation) + ' has not been imported or created in the project')

  # Set the default annoatations in the Mosaic project
  api_ps.setDefaultVariantAnnotations(mosaicConfig, projectId, annotationIds)

# Add a variant filter to a Mosaic project
def addVariantFilter(filepath, filename):

  # Get the information from the json file defining the filter
  try: filterFile = open(filepath + filename, 'r')
  except: fail('File ' + filepath + filename + ' could not be opened')
  try: data = json.load(filterFile)
  except: fail('File ' + filepath + filename + ' is not a valid json file')
  filterFile.close()

  # Get the name of the filter to create
  try: filterName = data['name']
  except: fail('Variant filter file, ' + str(filename) + ', defines a filter with no name. A name needs to be supplied')

# Set project attributes to indicate when and which version of Calypso has been run
def updateCalypsoAttributes(mosaicConfig, resourceVersion, projectAttributes, publicAttributes, version, projectId, api_pa):

  # Define the public attributes to import or update
  calypso = {}
  calypso['Calypso version']          = {'uid': False, 'id': False, 'value': str(version), 'inProject': False}
  calypso['Calypso date run']         = {'uid': False, 'id': False, 'value': str(date.today()), 'inProject': False}
  calypso['Calypso history']          = {'uid': False, 'id': False, 'value': str(version) + ':' + str(date.today()), 'inProject': False}
  calypso['Calypso resource version'] = {'uid': False, 'id': False, 'value': str(resourceVersion), 'inProject': False}

  # Loop over all the attributes in the project and look for the Calypso attributes
  for attributeUid in projectAttributes: 
    if projectAttributes[attributeUid]['name'] in calypso:
      attributeName  = projectAttributes[attributeUid]['name']
      attributeId    = projectAttributes[attributeUid]['id']
      attributeValue = projectAttributes[attributeUid]['value']
      calypso[attributeName]['uid']       = attributeUid
      calypso[attributeName]['id']        = attributeId
      calypso[attributeName]['inProject'] = True

      # If this is the Calypso history attribute, prepend the current value to the existing value, otherwise overwrite
      # the value
      if str(attributeName) == str('Calypso history'): calypso[attributeName]['value'] = str(calypso['Calypso history']['value']) + ',' + str(attributeValue)
      api_pa.updateProjectAttribute(mosaicConfig, projectId, calypso[attributeName]['id'], calypso[attributeName]['value'])
      calypso[attributeName]['inProject'] = True
  
  # Import all attributes that are not in the project
  for attribute in calypso:
    if not calypso[attribute]['inProject']:

      # Loop over the available public attributes and find the missing attribute
      isImported = False
      for publicAttributeUid in publicAttributes:
        if str(attribute) == publicAttributes[publicAttributeUid]['name']:
          attributeId = publicAttributes[publicAttributeUid]['id']
          api_pa.importProjectAttribute(mosaicConfig, projectId, attributeId, calypso[attribute]['value'])
          isImported = True
          break

      # If the attribute wasn't found, warn the user
      if not isImported: fail('Unable to import attribute ' + str(attribute) + '. Ensure this attribute exists in the public attribute project')

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
