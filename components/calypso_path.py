import os

# Build the path to allow importing additional modules
def buildPath(path, utils):

  # Add public-utils and required directories to the path
  path.append(utils)
  path.append(utils + 'common_components')
  path.append(utils + 'scripts')
  path.append(utils + 'api_commands')

  # return the updated paths
  return path

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
