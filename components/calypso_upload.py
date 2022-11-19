#!/usr/bin/python

from __future__ import print_function
from datetime import date

import json
import os

# Output a script to upload variants to Mosaic
def uploadVariants(workingDir, utils, configFile, projectId):

  # Open a script file
  uploadFileName = workingDir + 'calypso_upload_variants.sh'
  try: uploadFile = open(uploadFileName, 'w')
  except: fail('Could not open ' + str(uploadFileName) + ' to write to')

  # Write the command to file to upload the filtered variants
  print('# Upload variants to Mosaic', file = uploadFile)
  print('python ', utils, 'scripts/upload_variants.py \\', sep = '', file = uploadFile)
  print('  -c ', configFile + ' \\', sep = '', file = uploadFile)
  print('  -p ', str(projectId) + ' \\', sep = '', file = uploadFile)
  print('  -m "allele" \\', sep = '', file = uploadFile)
  print('  -i $FILTEREDVCF ', file = uploadFile)

  # Write the command to file to upload a variant set of the "rare disease" variants
  print('# Upload rare disease variants to Mosaic', file = uploadFile)
  print('python ', utils, 'scripts/upload_variants.py \\', sep = '', file = uploadFile)
  print('  -c ', configFile + ' \\', sep = '', file = uploadFile)
  print('  -p ', str(projectId) + ' \\', sep = '', file = uploadFile)
  print('  -m "allele" \\', sep = '', file = uploadFile)
  print('  -n "rare disease" \\', sep = '', file = uploadFile)
  print('  -i $RAREDISEASEVCF ', file = uploadFile)

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
