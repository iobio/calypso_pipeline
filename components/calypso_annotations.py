#!/usr/bin/python

from __future__ import print_function

import os

# Call the command to extract annotations from the filtered vcf into a tsv file
def createAnnotationTsv(mosaicInfo, resource, scriptDir, reference, configFile, mosaicJson, toolsDir, bashFile, files):
  tags = ''
  uids = ''
  outputFiles = []

  # Public resources will have a uid associated with the resource and so this can be extracted
  for annotation in mosaicInfo['resources'][resource]['annotations']:
    tags += ',' + annotation
    uids += ',' + mosaicInfo['resources'][resource]['annotations'][annotation]['uid']
  tags = tags.strip(',')
  uids = uids.strip(',')
  print(file = bashFile)
  print('  # Resource: ', resource, sep = '', file = bashFile)
  print('  echo -n "Creating tsv file for resource ', resource, '..."', sep = '', file = bashFile)

  # Generate a command for each independent VCF file created
  for i, annotateFile in enumerate(files):

    # Define the name of the output tsv file
    if len(files) == 1: outputFile = str(resource) + '.tsv'
    else: outputFile = str(resource) + '_' + str(i + 1) + '.tsv'
    outputFiles.append(outputFile)

    print('  python ', scriptDir, '/generate_annotation_tsv.py \\', sep = '', file = bashFile)
    print('    -c "', configFile, '" \\', sep = '', file = bashFile)
    print('    -e "', resource, '" \\', sep = '', file = bashFile)
    print('    -r "', reference, '" \\', sep = '', file = bashFile)
    print('    -m ', mosaicJson, ' \\', sep = '', file = bashFile)
    print('    -g "', tags, '" \\', sep = '', file = bashFile)
    print('    -d "', uids, '" \\', sep = '', file = bashFile)
    if toolsDir: print('    -s "', toolsDir, '" \\', sep = '', file = bashFile)
    print('    -i ', annotateFile, ' \\', sep = '', file = bashFile)
    print('    -o ', outputFile, sep = '', file = bashFile)
  print('  echo "complete"', file = bashFile)

  # Return the name of the output file
  return outputFiles

# Create the tsv files for SpliceAI annotations
def createSpliceAITsv(bashFile, scriptDir, toolsDir, configFile, mosaicResourseJson, reference, files):
  outputFiles = []

  print(file = bashFile)
  print('  # Resource: SpliceAI', sep = '', file = bashFile)
  print('  echo -n "Creating tsv file for resource SpliceAI..."', sep = '', file = bashFile)

  # Generate a command for each independent VCF file created
  for i, annotateFile in enumerate(files):

    # Define the name of the output file
    outputFile = 'spliceai_' + str(i + 1) + '.tsv'
    outputFiles.append(outputFile)

    print('  python ', scriptDir, '/generate_spliceai_tsv.py \\', sep = '', file = bashFile)
    print('    -c "', configFile, '" \\', sep = '', file = bashFile)
    print('    -m "', mosaicResourseJson, '" \\', sep = '', file = bashFile)
    print('    -r "', reference, '" \\', sep = '', file = bashFile)
    print('    -s "', toolsDir, '" \\', sep = '', file = bashFile)
    print('    -i ', annotateFile, ' \\', sep = '', file = bashFile)
    print('    -o ', outputFile, sep = '', file = bashFile)
  print('  echo "complete"', file = bashFile)

  # Return the name of the output file
  return outputFiles

# Call the command to extract HPO annotations from the filtered vcf into a tsv file
def createHpoTsv(hpoInfo, scriptDir, configFile, hpoString, projectId, genePhenotype, utils, toolsDir, bashFile, files):
  uids        = {}
  outputFiles = []

  # Get the uids for the HPO Overlaps and HPO Terms annotations
  for annotation in hpoInfo['annotations']: uids[annotation] = hpoInfo['annotations'][annotation]['uid']

  print(file = bashFile)
  print('  # Resource: HPO', sep = '', file = bashFile)
  print('  echo -n "Creating tsv file for resource HPO..."', sep = '', file = bashFile)

  # Generate a command for each independent VCF file created
  for i, annotateFile in enumerate(files):

    # Define the name of the output file
    outputFile = 'hpo_' + str(i + 1) + '.tsv'
    outputFiles.append(outputFile)

    print('  python ', scriptDir, '/hpo.py \\', sep = '', file = bashFile)
    print('    -c "', configFile, '" \\', sep = '', file = bashFile)
    print('    -r "', hpoString, '" \\', sep = '', file = bashFile)
    print('    -p ', projectId, ' \\', sep = '', file = bashFile)
    print('    -f "', genePhenotype, '" \\', sep = '', file = bashFile)
    print('    -l "', utils, '" \\', sep = '', file = bashFile)
    print('    -d "', uids['HPO Overlaps'], '" \\', sep = '', file = bashFile)
    print('    -e "', uids['HPO Terms'], '" \\', sep = '', file = bashFile)
    print('    -b "', uids['HPO Labels'], '" \\', sep = '', file = bashFile)
    if toolsDir: print('    -s "', toolsDir, '" \\', sep = '', file = bashFile)
    print('    -i ', annotateFile, ' \\', sep = '', file = bashFile)
    print('    -o ', outputFile, sep = '', file = bashFile)
  print('  echo "complete"', file = bashFile)

  # Return the name of the output file
  return outputFiles

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
