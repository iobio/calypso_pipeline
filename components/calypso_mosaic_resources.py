#!/usr/bin/python

from __future__ import print_function
import json
import os

# Read through the Mosaic json file describing the mosaic information for uploading annotations
def readMosaicJson(mosaicFilename, reference):
  mosaicInfo              = {}
  mosaicInfo['resources'] = {}

  # Try and open the file
  try: mosaicFile = open(mosaicFilename, 'r')
  except: fail('The file describing Mosaic related information (' + str(mosaicFilename) + ') could not be found')

  # Extract the json information
  try: mosaicData = json.loads(mosaicFile.read())
  except: fail('The json file (' + str(mosaicFilename) + ') is not valid')

  # Store the data version
  try: mosaicInfo['version'] = mosaicData['version']
  except: fail('The Mosaic json (' + str(mosaicFilename) + ') does not include a version')

  # Store the date the version was created
  try: mosaicInfo['date'] = mosaicData['date']
  except: fail('The Mosaic json (' + str(mosaicFilename) + ') does not include a date of creation')

  # Check that the resource json reference matches the selected reference
  try: mosaicReference = mosaicData['reference']
  except: fail('The Mosaic json does not include a reference genome')
  if reference != mosaicReference: fail('The selected reference (' + str(reference) + ') does not match the Mosaic json reference (' + str(mosaicReference) + ')')

  # Store information on annotationns to remove from the project, and defaults to set
  if 'remove' in mosaicData: mosaicInfo['remove'] = mosaicData['remove'] if 'remove' in mosaicData else []
  if 'default_annotations' in mosaicData: mosaicInfo['defaultAnnotations'] = mosaicData['default_annotations']

  # Loop over all the specified resources, and get all information required to pull the annotations
  # into Mosaic
  try: resources = mosaicData['resources']
  except: fail('The Mosaic json (' + str(mosaicFilename) + ') does not include a list of resource')
  for resource in resources:
    mosaicInfo['resources'][resource] = {}

    # The annotation must be identified as public or private. The other required information is dependent
    # on this
    try: mosaicInfo['resources'][resource]['annotation_type'] = resources[resource]['annotation_type']
    except: fail('The Mosaic json file does not include the "annotation_type" field for resource: ' + str(resource))

    # The resource "class" dictates how the annotated vcf will be processed
    try: mosaicInfo['resources'][resource]['class'] = resources[resource]['class']
    except: fail('The Mosaic json does not contain the "class" field for resource: ' + str(resource))

    # Some annotations need to be extracted from an INFO field that is not the name provided in the annotation. For example, the HGVS
    # codes need to be extracted from the CSQ field. The "info_field" provides information on where the info should be extracted from
    mosaicInfo['resources'][resource]['info_field'] = resources[resource]['info_field'] if 'info_field' in resources[resource] else False

    # Collect the information required only for public annotations
    if mosaicInfo['resources'][resource]['annotation_type'] == 'public':

      # Public annotations require the project id of the project that "owns" the annotation. Annotation values
      # need to be uploaded to the owner project
      try: project_id = resources[resource]['project_id']
      except: fail('The Mosaic json file does not include a "project_id" for public resource: ' + str(resource))
      mosaicInfo['resources'][resource]['project_id'] = project_id

    # Fail if the supplied annotation_type is neither public or private
    elif mosaicInfo['resources'][resource]['annotation_type'] != 'private': fail('Annotation_type for ' + str(resource) + ' must be "public" or "private"')

    # Loop over the annotation information
    mosaicInfo['resources'][resource]['annotations'] = {}
    for annotation in resources[resource]['annotations']:
      mosaicInfo['resources'][resource]['annotations'][annotation] = {}

      # Check that all required information is supplied and store it
      try: mosaicInfo['resources'][resource]['annotations'][annotation]['uid'] = resources[resource]['annotations'][annotation]['uid']
      except: fail('The Mosaic json does not contain the "uid" field for annotation "' + str(annotation) + '" for resource "' + str(resource) + '"')
      try: mosaicInfo['resources'][resource]['annotations'][annotation]['type'] = resources[resource]['annotations'][annotation]['type']
      except: fail('The Mosaic json does not contain the "type" field for annotation "' + str(annotation) + '" for resource "' + str(resource) + '"')

      # Check if the "is_sample" flag is set. This means there is an annotation for each project sample
      mosaicInfo['resources'][resource]['annotations'][annotation]['isSample'] = resources[resource]['annotations'][annotation]['is_sample'] if 'is_sample' in resources[resource]['annotations'][annotation] else False

      # Check if the "position" value is set. This defines the position in a compound annotation that the desired annotation can be found
      mosaicInfo['resources'][resource]['annotations'][annotation]['position'] = resources[resource]['annotations'][annotation]['position'] if 'position' in resources[resource]['annotations'][annotation] else False

  return mosaicInfo

# Loop over all of the public annotations defined in the Mosaic resources json file and check that they exist
# in Mosaic and can be imported
def checkPublicAnnotations(mosaicConfig, resources, publicAnnotations, projectId, api_va):
  for resource in resources:
    if resources[resource]['annotation_type'] == 'public':

      # Loop over all the annotations required for this resource
      for annotation in resources[resource]['annotations']:
        uid = resources[resource]['annotations'][annotation]['uid']

        # If the annotation uid is not in the publicAnnotations dictionary, then this annotation does not exist for
        # import. Fail and indicate that the resource file may contain errors
        if uid not in publicAnnotations: fail('Annotation with uid ' + str(uid) + ', for resource ' + str(resource) + ', is not available for import. Check that the uid in the Mosaic resources file is correct')

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
