from datetime import date
from os.path import exists
from sys import path

import argparse
import os
import sys
import calypso_resources as res

# Add the path of the common functions and import them
path.append(os.path.dirname(__file__) + '/components')
import tools_bcftools as bcftools
import calypso_mosaic_resources as mosr
import calypso_path as cpath

def main():
  global allowedClasses
  global bcftoolsExe

  # Parse the command line
  args = parseCommandLine()

  # Define the executable bcftools command
  if args.tools_directory:
    if not args.tools_directory.endswith('/'): args.tools_directory = str(args.tools_directory) + '/'
    bcftoolsExe = args.tools_directory + 'bcftools/bcftools'
  else: bcftoolsExe = 'bcftools'
  sys.path = cpath.buildPath(sys.path, args.utils_directory)
  import mosaic_config
  import api_variant_annotations as api_va

  # Parse the Mosaic config file to get the token and url for the api calls
  mosaicRequired = {'MOSAIC_TOKEN': {'value': args.token, 'desc': 'An access token', 'long': '--token', 'short': '-t'},
                    'MOSAIC_URL': {'value': args.url, 'desc': 'The api url', 'long': '--url', 'short': '-u'},
                    'MOSAIC_ATTRIBUTES_PROJECT_ID': {'value': args.attributes_project, 'desc': 'The public attribtes project id', 'long': '--attributes_project', 'short': '-a'}}
  mosaicConfig   = mosaic_config.mosaicConfigFile(args.config)
  mosaicConfig   = mosaic_config.commandLineArguments(mosaicConfig, mosaicRequired)

  # Read the mosaicJson file to get information on how to process different annotations
  mosaicInfo = mosr.readMosaicJson(args.mosaic_json, args.reference)

  # Open an output tsv file to write annotations to
  outputFile = open(args.output_tsv, 'w')

  # Get the annotations from Mosaic in order to get the uids. HGVS annotations are private and so wont' be available
  # in the resource file
  annotations       = api_va.getAnnotationDictNameIdUid(mosaicConfig, args.project_id)
  mosaicAnnotations = {}
  for annotationName in annotations: mosaicAnnotations[annotationName] = annotations[annotationName]['uid']

  # Loop over the annotations that are to be uploaded to Mosaic for this resource and get the annotation name, uid and
  # type
  annotations = {}
  uids        = []
  for annotation in mosaicInfo['resources']['HGVS']['annotations']:
    uid = mosaicAnnotations[annotation] if annotation in mosaicAnnotations else fail('Annotation ' + str(annotation) + ' is not present in the Mosaic project')
    tagType = mosaicInfo['resources']['HGVS']['annotations'][annotation]['type']
    annotations[annotation] = {'uid': uid, 'type': tagType}

  # Write the header line to the tsv file
  print('CHROM\tSTART\tEND\tREF\tALT\t', '\t'.join(str(x['uid']) for x in annotations.values()), sep = '', file = outputFile)

  # Check that the delimiter is provided to determine how to split up the compound annotation. If it is not set, then provide a
  # warning and do not preceed with this annotation
  try: delimiter = mosaicInfo['resources']['HGVS']['delimiter']
  except: fail('The delimiter field is not provided for resource ' + str(resource) + ' and so its annotation cannot be processed.')

  # Get all the MANE transcript ids
  maneFile = '/scratch/ucgd/lustre-work/marth/marth-projects/calypso/reannotation/data/GRCh38/reference/MANE.GRCh38.v1.3.ensembl.transcript_ids.txt'
  mane     = open(maneFile, 'r')
  maneIds  = []
  for record in mane.readlines(): maneIds.append(record.rstrip())
  mane.close()
  
  # Loop over all records in the vcf file
  for record in os.popen(bcftools.query(bcftoolsExe, args.input_vcf, [mosaicInfo['resources']['HGVS']['info_field']])).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')
  
    # If all values are '.', this line can be ignored
    uniqueValues = set(fields[5:])
    if len(uniqueValues) == 1 and list(uniqueValues)[0] == '.': continue

    # Update the chromosome and position
    fields[0], fields[2] = updateCoords(fields[0], fields[2])

    # If there are multiple annotations, e.g. there could be HGVS codes for multiple transcripts, break the
    # annotations up
    hgvsFields    = fields.pop()
    combinedCodes = hgvsFields.split(',') if ',' in hgvsFields else [hgvsFields]
    positions = []
    for annotation in mosaicInfo['resources']['HGVS']['annotations']:
      position = mosaicInfo['resources']['HGVS']['annotations'][annotation]['position']
      if position: positions.append(position)

    # Loop over the HGVS annotations
    options       = {}
    optionalCodes = {}
    for i, code in enumerate(combinedCodes):
      options[i]       = False
      hgvsCodes        = code.split(delimiter)
      optionalCodes[i] = {}
      for position in positions:
        fullCode = hgvsCodes[int(position) - 1]
        hgvsCode = fullCode.split(':')
        if not fullCode: continue
        if hgvsCode[0] in maneIds: options[i] = True
        optionalCodes[i][position] = hgvsCode
      if len(optionalCodes[i]) != 2 and options[i]: options[i] = 'mane'
      elif len(optionalCodes[i]) == 2 and not options[i]: options[i] = 'nonmane'
      elif len(optionalCodes[i]) == 0: del optionalCodes[i]

    # If no options were found (e.g. the HGVS codes were empty), skip this field
    if len(optionalCodes) == 0: continue

    # Loop over the available HGVS codes and choose the best option
    topChoices = []
    for i in options:
      if options[i] == True: topChoices.append(i)

    # Store all the HGVS codes as a string. The 'best' code will be provided by default, but all available
    # codes will be stored in a separate annotation
    allCDot = []
    allPDot = []
    for i in optionalCodes:
      if positions[0] in optionalCodes[i]:
        if optionalCodes[i][positions[0]][1] not in allCDot: allCDot.append(optionalCodes[i][positions[0]][1])
      if positions[1] in optionalCodes[i]:
        if optionalCodes[i][positions[1]][1] not in allPDot: allPDot.append(optionalCodes[i][positions[1]][1])

    # If there was a single top choice, use this:
    if len(topChoices) == 1:
      fields.append(optionalCodes[topChoices[0]][positions[0]][1])
      fields.append(optionalCodes[topChoices[0]][positions[1]][1])

    # If there are multiple top choices
    elif len(topChoices) > 1: fail('FAIL MULTIPLE TOP CHOICES')

    # If there are no top choices, look for MANE transcripts, but records that don't have both a c. and a p.
    else:
      maneChoices = []
      for i in options:
        if options[i] == 'mane': maneChoices.append(i)

      # If there was a single mane, use this, otherwise pick from the options if there are
      # multiple. If there are no manes, move to non-mane transcripts that have both a c. and p.
      if len(maneChoices) == 1:
        if positions[0] in optionalCodes[maneChoices[0]]: fields.append(optionalCodes[maneChoices[0]][positions[0]][1])
        else: fields.append('')
        if positions[1] in optionalCodes[maneChoices[0]]: fields.append(optionalCodes[maneChoices[0]][positions[1]][1])
        else: fields.append('')
      elif len(maneChoices) > 1:
        i = list(optionalCodes.keys())[0]
        if positions[0] in optionalCodes[i]: fields.append(optionalCodes[i][positions[0]][1])
        else: fields.append('')
        if positions[1] in optionalCodes[i]: fields.append(optionalCodes[i][positions[1]][1])
        else: fields.append('')
        #for position in positions:
        #  if position in optionalCodes[i]: fields.append(optionalCodes[i][position][1])
        #  else: fields.append('')
      else:
        nonManeChoices = []
        for i in options:
          if options[i] == 'nonmane': nonManeChoices.append(i)

        # Look through the options as before
        if len(nonManeChoices) == 1:
          if positions[0] in optionalCodes[nonManeChoices[0]]: fields.append(optionalCodes[nonManeChoices[0]][positions[0]][1])
          else: fields.append('')
          if positions[1] in optionalCodes[nonManeChoices[0]]: fields.append(optionalCodes[nonManeChoices[0]][positions[1]][1])
          else: fields.append('')
        elif len(nonManeChoices) > 1:
          i = list(optionalCodes.keys())[0]
          if positions[0] in optionalCodes[i]: fields.append(optionalCodes[i][positions[0]][1])
          else: fields.append('')
          if positions[1] in optionalCodes[i]: fields.append(optionalCodes[i][positions[1]][1])
          else: fields.append('')
        else:

          # If there is a single non mane, use it, otherwise just take the first value
          if len(optionalCodes) >= 1:
            i = list(optionalCodes.keys())[0]
            if positions[0] in optionalCodes[i]: fields.append(optionalCodes[i][positions[0]][1])
            else: fields.append('')
            if positions[1] in optionalCodes[i]: fields.append(optionalCodes[i][positions[1]][1])
            else: fields.append('')
  
    # Append the list of all available c. and p. values and write to file
    fields.append(', '.join(allCDot))
    fields.append(', '.join(allPDot))
    print('\t'.join(fields), file = outputFile)

  # Close the output tsv file
  outputFile.close()

# Input options
def parseCommandLine():
  global version
  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  parser.add_argument('--input_vcf', '-i', required = True, metavar = 'string', help = 'The input vcf file to annotate')
  parser.add_argument('--output_tsv', '-o', required = True, metavar = 'string', help = 'The output tsv file')
  parser.add_argument('--reference', '-r', required = True, metavar = 'string', help = 'The genome reference file used')
  parser.add_argument('--tools_directory', '-s', required = True, metavar = 'string', help = 'The path to the directory where the tools live')
  parser.add_argument('--utils_directory', '-l', required = True, metavar = 'string', help = 'The path to the public-utils directory')
  parser.add_argument('--mosaic_json', '-m', required = True, metavar = 'string', help = 'The json file describing the Mosaic parameters')
  parser.add_argument('--project_id', '-p', required = True, metavar = 'string', help = 'The project id that variants will be uploaded to')

  # Mosaic arguments
  parser.add_argument('--config', '-c', required = True, metavar = "string", help = "The config file for Mosaic")
  parser.add_argument('--token', '-t', required = False, metavar = "string", help = "The Mosaic authorization token")
  parser.add_argument('--url', '-u', required = False, metavar = "string", help = "The base url for Mosaic")
  parser.add_argument('--attributes_project', '-a', required = False, metavar = "integer", help = "The Mosaic project id that contains public attributes")

  # Version
  parser.add_argument('--version', '-v', action="version", version='Calypso annotation pipeline version: ' + str(version))

  return parser.parse_args()

# Update the chromosome and position in the tsv file
def updateCoords(chrom, pos):

  # Check that the chromosome does not include "chr" prefix
  if chrom.startswith('chr'): chrom = chrom.strip('chr')

  # Add one to the end position
  pos = str(int(pos) + 1)

  # Return the updated values
  return chrom, pos

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)

# Initialise global variables

# Pipeline version
version = "0.0.1"

# Mosaic aip information
mosaicConfig = {}

# Define the bcftools executable
bcftoolsExe = False

# Define the allowed annotation classes
allowedClasses = ['A', 'B', 'C', 'clinvar', 'compound', 'OMIM']

if __name__ == "__main__":
  main()