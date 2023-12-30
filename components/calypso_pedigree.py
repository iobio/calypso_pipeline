import os

# Determine if the project contains a proband, mother and father to allow detection of de novo
# and comp het variants
def isTrio(samples):
  hasMother = False
  hasFather = False
  for sample in samples:
    print(sample, samples[sample])
    if str(samples[sample]['relation']) == 'Mother':
      if hasMother: fail('The project has more than one sample defined as a mother')
      hasMother = True
    elif str(samples[sample]['relation']) == 'Father':
      if hasFather: fail('The project has more than one sample defined as a father')
      hasFather = True

  # Return true if the proband has a mother and a father with data
  if hasMother and hasFather: return True
  else: return False

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
