#!/usr/bin/python

from __future__ import print_function

import os

# Call the command to extract annotations from the filtered vcf into a tsv file
def createAnnotationTsv(mosaicInfo, resource, scriptDir, reference, configFile, mosaicJson, bashFile):
  tags = ''
  uids = ''

  # Public resources will have a uid associated with the resource and so this can be extracted
  for annotation in mosaicInfo['resources'][resource]['annotations']:
    tags += ',' + annotation
    uids += ',' + mosaicInfo['resources'][resource]['annotations'][annotation]['uid']
  tags = tags.strip(',')
  uids = uids.strip(',')

  # Define the name of the output tsv file
  outputFile = str(resource) + '.tsv'

  print(file = bashFile)
  print('  # Resource: ', resource, sep = '', file = bashFile)
  print('  echo -n "Creating tsv file for resource ', resource, '..."', sep = '', file = bashFile)
  print('  python ', scriptDir, '/generate_annotation_tsv.py \\', sep = '', file = bashFile)
  print('    -c "', configFile, '" \\', sep = '', file = bashFile)
  print('    -e "', resource, '" \\', sep = '', file = bashFile)
  print('    -r "', reference, '" \\', sep = '', file = bashFile)
  print('    -m ', mosaicJson, ' \\', sep = '', file = bashFile)
  print('    -g "', tags, '" \\', sep = '', file = bashFile)
  print('    -d "', uids, '" \\', sep = '', file = bashFile)
  print('    -i $FILTEREDVCF \\', sep = '', file = bashFile)
  print('    -o ', outputFile, sep = '', file = bashFile)
  print('  echo "complete"', file = bashFile)

  # Return the name of the output file
  return outputFile

# Call the command to extract HPO annotations from the filtered vcf into a tsv file
def createHpoTsv(hpoPublicInfo, hpoPrivateInfo, scriptDir, configFile, hpoString, projectId, genePhenotype, utils, bashFile):
  uids = {}

  # Get the uids for the HPO Overlaps and HPO Terms annotations
  for annotation in hpoPublicInfo['annotations']: uids[annotation] = hpoPublicInfo['annotations'][annotation]['uid']
  for annotation in hpoPrivateInfo['annotations']: uids[annotation] = hpoPrivateInfo['annotations'][annotation]['uid']

  print(file = bashFile)
  print('  # Resource: HPO', sep = '', file = bashFile)
  print('  echo -n "Creating tsv file for resource HPO..."', sep = '', file = bashFile)
  print('  python ', scriptDir, '/hpo.py \\', sep = '', file = bashFile)
  print('    -c "', configFile, '" \\', sep = '', file = bashFile)
  print('    -r "', hpoString, '" \\', sep = '', file = bashFile)
  print('    -p ', projectId, ' \\', sep = '', file = bashFile)
  print('    -o "', genePhenotype, '" \\', sep = '', file = bashFile)
  print('    -l "', utils, '" \\', sep = '', file = bashFile)
  print('    -d "', uids['HPO Overlaps'], '" \\', sep = '', file = bashFile)
  print('    -e "', uids['HPO Terms'], '" \\', sep = '', file = bashFile)
  print('    -b "', uids['HPO Labels'], '" \\', sep = '', file = bashFile)
  print('    -i $FILTEREDVCF', sep = '', file = bashFile)
  print('  echo "complete"', file = bashFile)

  # Return the name of the output file
  return 'hpo_public.tsv'

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
