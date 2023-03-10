from __future__ import print_function

import os
import json

# Determine the id of the proband
def getProband(mosaicConfig, args, workingDir, api_s, api_ped):
  samples         = {}
  mosaicSampleIds = {}
  deletePed       = False

  # Get the samples Mosaic id and store
  mosaicSamples = api_s.getSampleNamesAndIds(mosaicConfig, args.project_id)
  for sample in mosaicSamples: mosaicSampleIds[mosaicSamples[sample]] = sample

  # If a ped file was supplied on the command line, use this instead of the Mosaic information
  if args.ped: 

    # Open the ped file and get the samples
    try: pedFile = open(args.ped, 'r')
    except: fail('Couldn\'t open the ped file (' + args.ped + '). Components/calypso_samples.py (1).')

    # Loop over all lines in the ped file
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

        # Get the sex of the sample
        if int(fields[4]) == 1: sex = 'male'
        elif int(fields[4]) == 2: sex = 'female'
        else: sex = 'unknown'
  
        # Determine if the sample is affected
        try: affectedStatus = int(fields[5])
        except: fail('Ped file does not contain enough fields. Components/calypso_samples.py (3)')
        isAffected = True if int(fields[5]) == 2 else False
        samples[sample] = {'noParents': noParents, 'father': father, 'mother': mother, 'isAffected': isAffected, 'sex': sex}
  
        # Attach the Mosaic sample id
        try: samples[sample]['mosaicId'] = mosaicSamples[sample]
        except: fail('Sample ' + str(sample) + ' is not present in the Mosaic project. Components/calypso_samples.py (4)')
  
    # Close the ped file
    pedFile.close()

  # If no ped was supplied on the command, get the pedigree information from Mosaic
  else:
    for sample in mosaicSamples:
      for pedInfo in api_ped.getPedigree(mosaicConfig, args.project_id, mosaicSamples[sample]):
        sampleName = pedInfo['name']
        sampleId   = pedInfo['id']
  
        # Get information on the samples parents
        motherId   = pedInfo['pedigree']['maternal_id']
        fatherId   = pedInfo['pedigree']['paternal_id']
        if motherId and fatherId: noParents = 2
        elif not motherId and not fatherId: noParents = 0
        else: noParents = 1
  
        # Get the sample's sex
        if pedInfo['pedigree']['sex'] == 1: sex = 'male'
        elif pedInfo['pedigree']['sex'] == 2: sex = 'female'
        else: sex = 'unknown'
  
        # Determine if the sample is affected
        isAffected = True if pedInfo['pedigree']['affection_status'] == 2 else False
  
        # If this sample has already been seen, check that the information is as already recorded, otherwise
        # add to samples
        if sampleName in samples:
          if sampleId   != samples[sampleName]['mosaicId']: fail('Error in Mosaic pedigree information')
          if motherId   != samples[sampleName]['motherId']: fail('Error in Mosaic pedigree information')
          if fatherId   != samples[sampleName]['fatherId']: fail('Error in Mosaic pedigree information')
          if noParents  != samples[sampleName]['noParents']: fail('Error in Mosaic pedigree information')
          if sex        != samples[sampleName]['sex']: fail('Error in Mosaic pedigree information')
          if isAffected != samples[sampleName]['isAffected']: fail('Error in Mosaic pedigree information')
        else: samples[sampleName] = {'mosaicId': sampleId, 'motherId': motherId, 'fatherId': fatherId, 'noParents': noParents, 'sex': sex, 'isAffected': isAffected}

    # With all samples processed, populate the samples with the mother and father names as well as the ids
    for sample in samples:
      motherId = samples[sample]['motherId']
      fatherId = samples[sample]['fatherId']
      if motherId == None: samples[sample]['mother'] = None
      else: samples[sample]['mother'] = mosaicSampleIds[motherId]
      if fatherId == None: samples[sample]['father'] = None
      else: samples[sample]['father'] = mosaicSampleIds[fatherId]

  # With the pedigree information defined, determine which sample is the proband, if either parent is listed as
  # affected, and whether there are any affected / unaffected siblings

  # If no samples were found, fail.
  if len(samples) == 0: fail('No pedigree information was found for any samples')

  # Determine the proband and the family structure
  elif len(samples) == 1: samples, proband, familyType = singletonStructures(samples)
  elif len(samples) == 2: samples, proband, familyType = duoStructures(samples)
  elif len(samples) == 3: samples, proband, familyType = trioStructures(samples)
  #elif len(samples) == 4:

  # If there are more than 4 family members, additional logic needs to be added
  else: fail('Additional logic is required for families with more than 4 members')

  # If no ped was provided, write out a temporary ped file as the bash script will need one
  if not args.ped:
    deletePed = True
    args.ped  = str(workingDir) + 'calypso.ped'
    pedHandle = open(args.ped, 'w')
    print('#Kindred_ID\tSample_ID\tPaternal_ID\tMaternal_ID\tSex\tAffection_Status', file = pedHandle)
    for sample in mosaicSamples:
  
      # Define the sex and affectation_status
      if samples[sample]['sex'] == 'male': sex = 1
      elif samples[sample]['sex'] == 'female': sex = 2
      else: sex = 0
      if samples[sample]['isAffected']: affected = 2
      else: affected = 1
  
      # Get the names of the mother and father
      mother = samples[sample]['mother'] if samples[sample]['mother'] else 0
      father = samples[sample]['father'] if samples[sample]['father'] else 0
      print('Kindred', sample, father, mother, sex, affected, sep = '\t', file = pedHandle)

  # Return samples information
  return args, proband, samples, familyType

# Logic for determining structure of a singleton
def singletonStructures(samples):
  for sample in samples:

    # ...and they are unaffected, there is no proband, so fail
    if not samples[sample]['isAffected']: fail('Pedigree only contains a single, unaffected sample')

    # ...and parents are listed for the sample, the pedigree is incomplete, so also fail
    if samples[sample]['noParents'] != 0: fail('Pedigree contains a single sample, but parent ids are listed. The pedigree is therefore incomplete')

    # ...and they are affected, they are the proband.
    else:
      samples[sample]['relationship'] = 'Proband'
      proband                         = [samples[sample]['mosaicId']]
      familyType                      = 'singleton'

  # Return the updated samples along with the name of the proband and the family structure
  return samples, proband, familyType

# Logic for determining structure of a duo
def duoStructures(samples):

  # Loop over the 2 samples
  for i, sample in enumerate(samples):

    # Get information on the first sample
    if i == 0:
      name       = sample
      isAffected = samples[sample]['isAffected']
      mother     = samples[sample]['mother']
      father     = samples[sample]['father']
      noParents  = samples[sample]['noParents']

    # Get information on the second sample and determine the family structure. The following duo pedigrees are allowed:
    # 1. A proband and an unaffected sibling
    # 2. A proband and an affected sibling
    # 3. A proband and an unaffected parent
    # 4. A proband and an affected parent
    else:

      # If neither sample is affected, fail
      if not isAffected and not samples[sample]['isAffected']: fail('The pedigree contains 2 samples, but neither are listed as affected')

      # If both samples are affected
      elif isAffected and samples[sample]['isAffected']:

        # If the samples are siblings, they will have the same parents
        if mother == samples[sample]['mother'] and father == samples[sample]['father']:

          # Case 2. This pedigree is assumed to be two affected siblings. In this case, both samples are set as proband
          if samples[sample]['noParents'] == 0:
            samples[sample]['relationship'] = 'Proband'
            samples[name]['relationship']   = 'Proband'
            proband                         = [name, sample]
            familyType                      = 'two_affected_siblings'

          # If the samples are siblings, but the pedigree contains information about parents, the pedigree is incomplete as
          # there is no information about the parents
          else: fail('The pedigree contains 2 siblings with samples defined for parents, but no information about the parents is available')

        # For proband / parent duos, the parent cannot have any parents of their own defined
        else:

          # If both samples have parents listed, fail as multi generational families are not handled
          if noParents != 0 and samples[sample]['noParents'] != 0: fail('Pedigree is a duo, but both samples have parents listed')

          # The 2 samples don't share the both same parents, so if they share either they are step-siblings which are not
          # handled
          if mother == samples[sample]['mother'] or father == samples[sample]['father']: fail('Pedigree is a duo of step-siblings which is not handled')

          # Case 4 could include an affected mother or father
          if name == samples[sample]['mother']:
            samples[sample]['relationship'] = 'Proband'
            samples[name]['relationship']   = 'Mother'
            proband                         = [sample]
            familyType                      = 'duo_affected_mother'
          elif name == samples[sample]['father']:
            samples[sample]['relationship'] = 'Proband'
            samples[name]['relationship']   = 'Father'
            proband                         = [sample]
            familyType                      = 'duo_affected_father'
          elif sample == mother:
            samples[sample]['relationship'] = 'Mother'
            samples[name]['relationship']   = 'Proband'
            proband                         = [name]
            familyType                      = 'duo_affected_mother'
          elif sample == father:
            samples[sample]['relationship'] = 'Father'
            samples[name]['relationship']   = 'Proband'
            proband                         = [name]
            familyType                      = 'duo_affected_father'
          else: fail('Unknown duo structure')

      # If a single sample is affected, this must be case 1 or 3
      else:

        # If the two samples share the same parents, they are either siblings, or, if neither has any parents listed,
        # they could be 2 unrelated individuals. However, assume these are siblings in this case
        if mother == samples[sample]['mother'] and father == samples[sample]['father']:

          # Case 1. Siblings
          if noParents == 0:
            if isAffected: 
              samples[sample]['relationship'] = 'Sibling'
              samples[name]['relationship']   = 'Proband'
              proband                         = [name]
              familyType                      = 'duo_sibling'
            else:
              samples[sample]['relationship'] = 'Proband'
              samples[name]['relationship']   = 'Sibling'
              proband                         = [sample]
              familyType                      = 'duo_sibling'

          # Case 1
          else:
            if isAffected: 
              samples[sample]['relationship'] = 'Sibling'
              samples[name]['relationship']   = 'Proband'
              proband                         = [name]
              familyType                      = 'duo_sibling'
            else:
              samples[sample]['relationship'] = 'Proband'
              samples[name]['relationship']   = 'Sibling'
              proband                         = [sample]
              familyType                      = 'duo_sibling'

        # If the samples do not have the same parents, this must be case 3
        else:
          if name == samples[sample]['mother']:
            if isAffected: fail('The defined pedigree (affected mother and unaffected child) is not supported')
            samples[sample]['relationship'] = 'Proband'
            samples[name]['relationship']   = 'Mother'
            proband                         = [sample]
            familyType                      = 'duo_unaffected_mother'
          elif name == samples[sample]['father']:
            if isAffected: fail('The defined pedigree (affected father and unaffected child) is not supported')
            samples[sample]['relationship'] = 'Proband'
            samples[name]['relationship']   = 'Father'
            proband                         = [sample]
            familyType                      = 'duo_unaffected_father'
          elif sample == mother:
            if not isAffected: fail('The defined pedigree (affected mother and unaffected child) is not supported')
            samples[sample]['relationship'] = 'Mother'
            samples[name]['relationship']   = 'Proband'
            proband                         = [name]
            familyType                      = 'duo_unaffected_mother'
          elif sample == father:
            if not isAffected: fail('The defined pedigree (affected father and unaffected child) is not supported')
            samples[sample]['relationship'] = 'Father'
            samples[name]['relationship']   = 'Proband'
            proband                         = [name]
            familyType                      = 'duo_unaffected_father'
          else: fail('Unknown duo structure')
  
  # Return the updated samples along with the name of the proband and the family structure
  return samples, proband, familyType

# Logic for determining structure of a trio
def trioStructures(samples):

  # Loop over the 3 samples and determine the number of affected samples, and the number of samples with parents
  noAffecteds   = 0
  noWithParents = 0
  for sample in samples:
    if samples[sample]['isAffected']: noAffecteds += 1
    if samples[sample]['noParents'] != 0: noWithParents += 1

  # If none of the samples have parents, this is assumed to correspond to three siblings. Determine which sibs are
  # affected and set the proband accordingly
  if noWithParents == 0:
    proband = []
    for sample in samples:
      if samples[sample]['isAffected']:
        proband.append(sample)
        samples[sample]['relationship'] = 'Proband'
      else: samples[sample] = 'sibling'
    if len(proband) == 0: fail('The pedigree consists of multiple unaffected siblings. At least one sample must be affected')
    elif len(proband) == 1: familyType = 'one_affected_two_unaffected_siblings'
    elif len(proband) == 2: familyType = 'two_affected_one_unaffected_siblings'
    elif len(proband) == 3: familyType = 'three_affected_siblings'

  # If there is a single sample with parents, check that the parents are correct. If so, this is a single
  # proband with unaffected parents
  elif noWithParents == 1:

    # Loop over the samples, find the sample with parents and check they are the other two samples in the pedigree
    for sample in samples:
      if samples[sample]['noParents'] == 0: continue
      else:
        mother = samples[sample]['mother']
        father = samples[sample]['father']
        if mother not in samples: fail('Proband in trio does not have the mother listed in the pedigree')
        if father not in samples: fail('Proband in trio does not have the father listed in the pedigree')

        # Ensure the parental sex's are correct
        if samples[mother]['sex'] != 'female': fail('The mother of the proband in a family trio is not listed as female in the pedigree')
        if samples[father]['sex'] != 'male': fail('The father of the proband in a family trio is not listed as male in the pedigree')

        # This sample is the proband, so must be listed as affected
        if not samples[sample]['isAffected']: fail('The child in a family trio is not listed as affected, which is not handled')

        # Update the family information
        samples[sample]['relationship'] = 'Proband'
        samples[mother]['relationship'] = 'Mother'
        samples[father]['relationship'] = 'Father'
        proband                         = [sample]
        familyType                      = 'trio'

  # Not handled other structures of parents
  else: fail('Unknown trio structure - logic needs updating to accommodate')

  # Return the updated samples along with the name of the proband and the family structure
  return samples, proband, familyType

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
