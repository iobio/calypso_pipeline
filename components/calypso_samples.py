#!/usr/bin/python

from __future__ import print_function

import os
import json

# Determine the id of the proband
def getProband(mosaicConfig, ped, familyType, projectId, api_s):
  samples = {}

  # Open the ped file and get the samples
  try: pedFile = open(ped, "r")
  except: fail('Couldn\'t open the ped file (' + ped + '). Components/calypso_samples.py (1).')

  # Get the samples Mosaic id and store
  mosaicSamples = {}
  try: data = json.loads(os.popen(api_s.getSamples(mosaicConfig, projectId)).read())
  except: fail('Failed to read samples (components/calypso_samples.py)')

  # If the API call returned an error, fail and print the error message
  if 'message' in data: fail('Failure in components/calypso_samples.py (2). Message: ' + data['message'])
  for sample in data: mosaicSamples[sample['name']] = sample['id']

  # Get information on all the samples
  noAffected = 0
  for line in pedFile:

    # Ignore the header line
    if not line.startswith('#'):
      fields = line.rstrip().split('\t')
      sample = fields[1]

      # Determine if the sample has parents
      father = False if fields[2] == '0' else fields[2]
      mother = False if fields[3] == '0' else fields[3]
      if mother and father: noParents = 2
      elif not mother and not father: noParents = 0
      else: noParents = 1

      # Determine if the sample is affected
      try: affectedStatus = int(fields[5])
      except: fail('Ped file does not contain enough fields. Components/calypso_samples.py (3)')
      if int(fields[5]) == 2:
        isAffected = True 
        noAffected += 1
        proband = sample
      else: isAffected = False
      samples[sample] = {'noParents': noParents, 'father': father, 'mother': mother, 'isAffected': isAffected, 'relationship': False}
      if isAffected: samples[sample]['relationship'] = 'Proband'

      # Attach the Mosaic sample id
      try: samples[sample]['mosaicId'] = mosaicSamples[sample]
      except: fail('Sample ' + str(sample) + ' is not present in the Mosaic project. Components/calypso_samples.py (4)')

  # Close the ped file
  pedFile.close()

  # Check that the ped file conforms to the supplied family type
  if familyType == 'singleton':
    if len(samples) != 1: fail('Family type was specified as a singleton, but the ped file contains multiple samples')
  elif familyType == 'duo':
    if len(samples) != 2: fail('Family type was specified as a duo, but the ped file doesn\'t contain two samples')
  elif familyType == 'trio':
    if len(samples) != 3: fail('Family type was specified as a trio, but the ped file doesn\'t contain three samples')
  elif familyType == 'quad':
    if len(samples) != 4: fail('Family type was specified as a quad, but the ped file doesn\'t contain four samples')

  # If multiple samples are affected, throw a warning
  if noAffected != 1: fail('Cannot determine proband - multiple samples in the ped are affected')

  # Identify the mother and father of the proband
  mother = samples[proband]['mother']
  father = samples[proband]['father']
  if samples[proband]['mother']: samples[mother]['relationship'] = 'Mother'
  if samples[proband]['father']: samples[father]['relationship'] = 'Father'

  # Identify siblings
  for sample in samples:
    if samples[sample]['mother'] == mother and samples[sample]['father'] and not samples[sample]['isAffected']: samples[sample]['relationship'] = 'Sibling'
    if not samples[sample]['relationship']: fail('Sample ' + str(sample) + ' has an unknown relationship to the proband')

  # Return samples information
  return proband, samples

# Get the order that the samples appear in, in the vcf header
def getSampleOrder(samples, vcf):

  # Get the final head line
  data = os.popen('bcftools view -h ' + str(vcf) + ' | tail -1').read().rstrip().split("\t")

  # Read through the samples and define the sample order
  for index, sample in enumerate(data[9:]):
    if sample in samples: samples[sample]['vcf_position'] = index

  # Check that every sample in samples has vcf_position set
  for sample in samples:
    if 'vcf_position' not in samples[sample]: fail('Sample ' + str(sample) + 'is listed in the ped file, but does not appear in the vcf header')

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
