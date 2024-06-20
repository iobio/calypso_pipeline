from datetime import date
from os.path import exists
from sys import path
from pprint import pprint

import argparse
import os
import sys

# Add the path of the common functions and import them
path.append(os.path.dirname(__file__) + '/components')
import tools_bcftools as bcftools
import calypso_mosaic_resources as mosr
import calypso_resources as res

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

  # Read the mosaicJson file to get information on how to process different annotations
  mosaicInfo = mosr.readMosaicJson(args.mosaic_json, args.reference)

  # Open an output tsv file to write annotations to
  outputFile = open('hgvs.tsv', 'w')

  # Loop over the annotations that are to be uploaded to Mosaic for this resource and get the annotation name, uid and
  # type
  annotations = {}
  uids        = []
  for annotation in mosaicInfo['resources']['HGVS']['annotations']:
    uid     = mosaicInfo['resources']['HGVS']['annotations'][annotation]['uid']
    tagType = mosaicInfo['resources']['HGVS']['annotations'][annotation]['type']
    annotations[annotation] = {'uid': uid, 'type': tagType}

  # Write the header line to the tsv file
  print('CHROM\tSTART\tEND\tREF\tALT\t', '\t'.join(str(x['uid']) for x in annotations.values()), sep = '', file = outputFile)

  # Check that the delimiter is provided to determine how to split up the compound annotation. If it is not set, then provide a
  # warning and do not preceed with this annotation
  try:
    delimiter = mosaicInfo['resources']['HGVS']['delimiter']
  except:
    fail('The delimiter field is not provided for resource ' + str(resource) + ' and so its annotation cannot be processed.')

  # Get all the MANE transcript ids
  maneFile = '/scratch/ucgd/lustre-work/marth/marth-projects/calypso/reannotation/data/GRCh38/reference/MANE.GRCh38.v1.3.ensembl.transcript_ids.txt'
  mane     = open(maneFile, 'r')
  maneIds  = []
  for record in mane.readlines():
    maneIds.append(record.rstrip())
  mane.close()
  
  # Loop over all records in the vcf file
  for record in os.popen(bcftools.query(bcftoolsExe, args.input_vcf, [mosaicInfo['resources']['HGVS']['info_field']])).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')
  
    # If all values are '.', this line can be ignored
    uniqueValues = set(fields[5:])
    if len(uniqueValues) == 1 and list(uniqueValues)[0] == '.':
      continue

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
        if '%' in hgvsCode[1]: hgvsCode[1] = hgvsCode[1].replace('%3D', '=')
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

    # If there was a single top choice, use this:
    if len(topChoices) == 1:
      fields = updateFields(fields, optionalCodes[topChoices[0]][positions[0]][1])
      fields = updateFields(fields, optionalCodes[topChoices[0]][positions[1]][1])

    # If there are multiple top choices, include all of them. These will likely come from different, overlapping genes
    elif len(topChoices) > 1:
      code = ''
      for i in topChoices: code += str(optionalCodes[i][positions[0]][1]) + ','
      code.rstrip(',')
      fields = updateFields(fields, code)

      code = ''
      for i in topChoices: code += str(optionalCodes[i][positions[1]][1]) + ','
      code.rstrip(',')
      fields = updateFields(fields, code)

    # If there are no top choices, look for MANE transcripts, but records that don't have both a c. and a p.
    else:
      maneChoices = []
      for i in options:
        if options[i] == 'mane': maneChoices.append(i)

      # If there was a single mane, use this, otherwise pick from the options if there are
      # multiple. If there are no manes, move to non-mane transcripts that have both a c. and p.
      if len(maneChoices) == 1:
        if positions[0] in optionalCodes[maneChoices[0]]: fields = updateFields(fields, optionalCodes[maneChoices[0]][positions[0]][1])
        else: fields.append('')
        if positions[1] in optionalCodes[maneChoices[0]]: fields = updateFields(fields, optionalCodes[maneChoices[0]][positions[1]][1])
        else: fields.append('')
      elif len(maneChoices) > 1:
        i = list(optionalCodes.keys())[0]
        if positions[0] in optionalCodes[i]: fields = updateFields(fields, optionalCodes[i][positions[0]][1])
        else: fields.append('')
        if positions[1] in optionalCodes[i]: fields = updateFields(fields, optionalCodes[i][positions[1]][1])
        else: fields.append('')
      else:
        nonManeChoices = []
        for i in options:
          if options[i] == 'nonmane': nonManeChoices.append(i)

        # Look through the options as before
        if len(nonManeChoices) == 1:
          if positions[0] in optionalCodes[nonManeChoices[0]]: fields = updateFields(fields, optionalCodes[nonManeChoices[0]][positions[0]][1])
          else: fields.append('')
          if positions[1] in optionalCodes[nonManeChoices[0]]: fields = updateFields(fields, optionalCodes[nonManeChoices[0]][positions[1]][1])
          else: fields.append('')
        elif len(nonManeChoices) > 1:
          i = list(optionalCodes.keys())[0]
          if positions[0] in optionalCodes[i]: fields = updateFields(fields, optionalCodes[i][positions[0]][1])
          else: fields.append('')
          if positions[1] in optionalCodes[i]: fields = updateFields(fields, optionalCodes[i][positions[1]][1])
          else: fields.append('')
        else:

          # If there is a single non mane, use it, otherwise just take the first value
          if len(optionalCodes) >= 1:
            i = list(optionalCodes.keys())[0]
            if positions[0] in optionalCodes[i]: fields = updateFields(fields, optionalCodes[i][positions[0]][1])
            else: fields.append('')
            if positions[1] in optionalCodes[i]: fields = updateFields(fields, optionalCodes[i][positions[1]][1])
            else: fields.append('')
  
    print('\t'.join(fields), file = outputFile)

  # Close the output tsv file
  outputFile.close()

# Input options
def parseCommandLine():
  global version
  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  parser.add_argument('--input_vcf', '-i', required = True, metavar = 'string', help = 'The input vcf file to annotate')
  parser.add_argument('--reference', '-r', required = True, metavar = 'string', help = 'The genome reference file used')
  parser.add_argument('--tools_directory', '-s', required = True, metavar = 'string', help = 'The path to the directory where the tools live')
  parser.add_argument('--mosaic_json', '-m', required = True, metavar = 'string', help = 'The json file describing the Mosaic parameters')

  return parser.parse_args()

# Update the chromosome and position in the tsv file
def updateCoords(chrom, pos):

  # Check that the chromosome does not include "chr" prefix
  if chrom.startswith('chr'): chrom = chrom.strip('chr')

  # Add one to the end position
  pos = str(int(pos) + 1)

  # Return the updated values
  return chrom, pos

# Add the value to the fields list ensuring that it is a valid value
def updateFields(fields, value):

  # Ensuret he value is under the 255 character limit
  if len(value) > 254: value = 'HGVS code too long'

  # Append the value and return
  fields.append(value)
  return fields

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)

# Initialise global variables

# Define the bcftools executable
bcftoolsExe = False

if __name__ == "__main__":
  main()
