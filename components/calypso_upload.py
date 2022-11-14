#!/usr/bin/python

from __future__ import print_function
from datetime import date

import json
import os

# Output a script to upload variants to Mosaic
def uploadVariants(mosaicInfo, mosaicConfig, workingDir, config, vcf, projectId, api_v):

  # Open a script file
  uploadFileName = workingDir + 'calypso_upload_variants.sh'
  try: uploadFile = open(uploadFileName, 'w')
  except: fail('Could not open ' + str(uploadFileName) + ' to write to')

  # Define the name of the variant set that will be create
  #description = 'Calypso_v' + str(mosaicInfo['version']) + '_variants_' + str(date.today())

  # The command to be executed requires the Mosaic token which should not be written in the script. Since the
  # config file has been validated, it contains a line of the form MOSAIC_TOKEN=X. Extract the token into a 
  # variable that can be used
  print('TOKEN=`less ', config, ' | grep MOSAIC_TOKEN`', sep = '', file = uploadFile)
  print(file = uploadFile)

  # Write the command to file
  print('# Upload variants to Mosaic', file = uploadFile)

  # Get the api command to upload variants
  mosaicConfig['token'] = '$TOKEN'
  command = api_v.postUploadVariants(mosaicConfig, vcf, "allele", projectId)
  print(command, file = uploadFile)

  # Close the file
  uploadFile.close()

  # Make the annotation script executable
  makeExecutable = os.popen('chmod +x ' + str(uploadFileName)).read()

# Open a script file
def openUploadAnnotationsFile(workingDir):
  uploadFilename = workingDir + 'calypso_upload_annotations.sh'
  try: uploadFile = open(uploadFilename, 'w')
  except: fail('Could not open ' + str(uploadFilename) + ' to write to')

  return uploadFilename, uploadFile

# Write out the command to upload an annotations tsv to Mosaic
def uploadAnnotations(utils, tsvFile, projectId, configFile, uploadFile):

  # Print the command to upload annotations to the output file
  print(file = uploadFile)
  print('# Upload ', tsvFile, ' annotations', sep = '', file = uploadFile)
  print('python ', utils, 'scripts/upload_annotations.py \\', sep = '', file = uploadFile)
  print('  -i ', tsvFile, ' \\', sep = '', file = uploadFile)
  print('  -p ', projectId, ' \\', sep = '', file = uploadFile)
  print('  -c ', configFile, sep = '', file = uploadFile)

# Close the upload annotations file
def closeUploadAnnotationsFile(uploadFilename, uploadFile):

  # Close the file
  uploadFile.close()

  # Make the annotation script executable
  makeExecutable = os.popen('chmod +x ' + str(uploadFilename)).read()

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
