#!/usr/bin/python

from __future__ import print_function
from datetime import date

import json
import os

# Output a script to upload variants to Mosaic
def uploadVariants(mosaicInfo, mosaicConfig, workingDir, config, vcf, projectId, api_v):

  # Open a script file
  uploadFileName = workingDir + 'calypso_upload_variants_to_mosaic.sh'
  try: uploadFile = open(uploadFileName, 'w')
  except: fail('Could not open ' + str(uploadFileName) + ' to write to')

  # Define the name of the variant set that will be create
  description = 'Calypso_v' + str(mosaicInfo['version']) + '_variants_' + str(date.today())

  # The command to be executed requires the Mosaic token which should not be written in the script. Since the
  # config file has been validated, it contains a line of the form MOSAIC_TOKEN=X. Extract the token into a 
  # variable that can be used
  print('TOKEN=`less ', config, ' | grep MOSAIC_TOKEN`', sep = '', file = uploadFile)
  print(file = uploadFile)

  # Write the command to file
  print('# Upload variants to Mosaic', file = uploadFile)

  # Get the api command to upload variants
  mosaicConfig['token'] = '$TOKEN'
  command = api_v.postUploadVariants(mosaicConfig, vcf, "allele", description, projectId)
  print(command, file = uploadFile)

  # Close the file
  uploadFile.close()

  # Make the annotation script executable
  makeExecutable = os.popen('chmod +x ' + str(uploadFileName)).read()

# Output scripts to upload annotations to Mosaic
#def uploadAnnotations():
#  global tsvFiles
#  global mosaicInfo
#  global workingDir
#  uploadFileName  = workingDir + "calypso_upload_annotations_to_mosaic.sh"
#  isScript        = False
#
#  # Loop over all resources
#  for resource in tsvFiles:
#    annotationType = mosaicInfo["resources"][resource]["annotation_type"]
#    isComplete     = False
#
#    # Public annotations need to be imported into the Mosaic project
#    uids = {}
#    if annotationType == "public": isComplete = getPublicAnnotation(args, resource)
#    elif annotationType == "private": isComplete = True
#
#    # Only create the script file if the annotation(s) exists in Mosaic
#    if isComplete: 
#
#      # Create a single script file to upload all variants
#      if not isScript: 
#        try: uploadFile = open(uploadFileName, "w")
#        except: fail("Could not open " + str(uploadFileName) + " to write to")
#        isScript = True
#
#      # Write the command to file
#      print("# Upload ", resource, " annotations to Mosaic", file = uploadFile)
#      print('echo -n "Uploading ', resource, '..."', sep = "", file = uploadFile)
#      print("python ", os.path.dirname( __file__ ), "/mosaic_commands/upload_annotations_to_mosaic.py \\", sep = "", file = uploadFile)
#      print("  -i ", str(tsvFiles[resource]), " \\", sep = "", file = uploadFile)
#      if args.config: print("  -c", str(args.config), "\\", file = uploadFile)
#      else: print("  -c \"Insert config file here\" \\", sep = "", file = uploadFile)
#  
#      # Public annotations need to be uploaded to the project that contains them, private
#      # annotations are uploaded to the defined project.
#      if annotationType == "private":
#        if args.project_id: print("  -p", str(args.project_id), file = uploadFile)
#        else: print("  -p \"Insert Mosaic project id here\"", file = uploadFile)
#      else: print("  -p ", mosaicInfo["resources"][resource]["project_id"], sep = "", file = uploadFile)
#      print('echo "complete"', file = uploadFile)
#      print(file = uploadFile)
#  
#  # If a script files was created, close it and make it executable
#  if isScript:
#
#    # Close the file
#    uploadFile.close()
#  
#    # Make the annotation script executable
#    makeExecutable = os.popen("chmod +x " + str(uploadFileName)).read()

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
