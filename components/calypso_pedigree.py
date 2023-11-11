import os

# Generate a ped file
def generatePedFile(api_ped, mosaicConfig, args, workingDir, proband, sampleId):

  # Open a ped file in the working directory
  pedFileName = str(workingDir) + str(proband) + '.ped'
  pedFile     = open(pedFileName, 'w')

  # Get the pedigree information, write to file and set args.ped to this file
  for line in api_ped.getPedLines(mosaicConfig, args.project_id, sampleId): print(line, file = pedFile)
  args.ped = pedFileName

  # Close the ped file
  pedFile.close()

  # Return the updated args
  return args

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
