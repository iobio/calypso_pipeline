from datetime import date

import json
import os

# Output a script to upload variants to Mosaic
def uploadVariants(workingDir, utils, configFile, projectId, vcfFiles):

  # Open a script file
  uploadFileName = workingDir + 'calypso_upload_variants.sh'
  try: uploadFile = open(uploadFileName, 'w')
  except: fail('Could not open ' + str(uploadFileName) + ' to write to')

  # Loop over all the created vcf files and write the command to upload them
  for vcf in vcfFiles:

    # Write the command to file to upload the filtered variants
    print('# Upload variants to Mosaic', file = uploadFile)
    print('API_CLIENT=/uufs/chpc.utah.edu/common/HIPAA/u0991467/programs/api_client', file = uploadFile)
    print('CONFIG=/uufs/chpc.utah.edu/common/HIPAA/u0991467/programs/configs/config_utah.ini', file = uploadFile)
    print('python3 $API_CLIENT/variants/upload_variants.py \\', sep = '', file = uploadFile)
    print('  -a $API_CLIENT \\', sep = '', file = uploadFile)
    print('  -c $CONFIG \\', sep = '', file = uploadFile)
    print('  -p ', str(projectId) + ' \\', sep = '', file = uploadFile)
    print('  -m "no-validation" \\', sep = '', file = uploadFile)
    print('  -v ', vcf, sep = '', file = uploadFile)
    print(file = uploadFile)

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
  print('echo -n "Uploading annotations from file ', str(tsvFile), '..."', sep = '', file = uploadFile)
  print('python3 ', utils, 'scripts/upload_annotations.py \\', sep = '', file = uploadFile)
  print('  -i ', tsvFile, ' \\', sep = '', file = uploadFile)
  print('  -p ', projectId, ' \\', sep = '', file = uploadFile)
  print('  -c ', configFile, sep = '', file = uploadFile)
  print('echo "complete"', file = uploadFile)

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
