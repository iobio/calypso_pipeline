from os.path import exists
from sys import path

import argparse
import os
import sys

# Add the path of the common functions and import them
scriptFile = os.path.dirname(__file__)
path.append(os.path.dirname(scriptFile) + '/components')
import tools_bcftools as bcftools
import calypso_mosaic_resources as mosr

def main():

  # Parse the command line
  args = parseCommandLine()

  # Define the executable bcftools command
  if args.tools_directory:
    if not args.tools_directory.endswith('/'): args.tools_directory = str(args.tools_directory) + '/'
    bcftoolsExe = args.tools_directory + 'bcftools/bcftools'
  else: bcftoolsExe = 'bcftools'

  # Check that the exectuable exists
  isExecutable = bcftools.isExecutable(bcftoolsExe)
  if not isExecutable: fail('bcftools executable (' + str(bcftoolsExe) + ') is not executable. Use --tools_directory (-s) to define the directory where the bcftools executable can be found')

  # Read the mosaicJson file to get information on how to process different annotations
  mosaicInfo = mosr.readMosaicJson(args.mosaic_json, args.reference)

  # Get the SpliceAI information from the Mosaic resource json
  if 'resources' not in mosaicInfo: fail('No "resources" section in the Mosaic resources json')
  if 'SpliceAI' not in mosaicInfo['resources']: fail('No "SpliceAI" information in the "resources" section of the Mosaic resources json')

  # Loop over the SpliceAI annotations to include and generate an ordered list of SpliceAI annotations
  annotationList      = []
  annotationUids      = ''
  annotationPositions = {}
  for annotation in mosaicInfo['resources']['SpliceAI']['annotations']:
    annotationName = mosaicInfo['resources']['SpliceAI']['annotations'][annotation]['fields']

    # Exclude the max score from the list
    if str(annotation) != 'SpliceAI Max Score':
      annotationList.append(annotationName)
      annotationUids += str(mosaicInfo['resources']['SpliceAI']['annotations'][annotation]['uid']) + '\t'

    # Store the uid of the max score to add to the end of the list of uids
    else:
      maxScoreUid = mosaicInfo['resources']['SpliceAI']['annotations'][annotation]['uid']

  # Add the max score uid to the end of the uid list
  annotationUids += maxScoreUid

  # Open an output tsv file to write annotations to
  outputFile = open(args.output_tsv, 'w')
  
  # Write the header line to the tsv file
  print('CHROM\tSTART\tEND\tREF\tALT\t', annotationUids, sep = '', file = outputFile)

  # Loop over all records in the vcf file
  for record in os.popen(bcftools.query(bcftoolsExe, args.input_vcf, annotationList)).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')

    # If all values are '.', this line can be ignored
    uniqueValues = set(fields[5:])
    if len(uniqueValues) == 1 and list(uniqueValues)[0] == '.': continue

    # Update the chromosome and position. If the chromosome begins with "chr", this needs to be stripped off. The end coordinate
    # (fields[2] - the start coordinate is fields[1]) needs to have one added to create a 0-based, half open coordinate
    fields[0], fields[2] = updateCoords(fields[0], fields[2])

    # The output line for each vcf record is the first 5 fields (coordinates and alleles), then the annotations in order.
    # Create the start of this line
    output = '\t'.join(fields[0:5])

    # The SpliceAI annotations appear in the order set in annotationList (as this was used to generate the bcftools query
    # command). Loop over the annotations in this order, if there are multiple values, use the maximum value and determine
    # the maximum score value across all the annotations
    i        = 5
    maxScore = 0
    for annotation in annotationList:
      maxValue = 0
      if ',' in fields[i]:
        for value in fields[i].split(','):
          if abs(float(value)) > abs(maxValue): maxValue = float(value)
      else: maxValue = fields[i]
      output += '\t' + str(maxValue)

      # Check if this is the max value seen so far of the SpliceAI scores
      if annotation.startswith('SpliceAI_DS'):
        if abs(float(maxValue)) > abs(maxScore): maxScore = float(maxValue)
      i += 1

    # Include the max score as the last annotation
    output += '\t' + str(maxScore)

    # Write all the annotations to the output tsv
    print(output, file = outputFile)
  
  # Close the output tsv file
  outputFile.close()

# Input options
def parseCommandLine():
  global version
  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  parser.add_argument('--config', '-c', required = True, metavar = "string", help = "The config file for Mosaic")
  parser.add_argument('--input_vcf', '-i', required = True, metavar = 'string', help = 'The input vcf file to annotate')
  parser.add_argument('--output_tsv', '-o', required = True, metavar = 'string', help = 'The output tsv file')
  parser.add_argument('--reference', '-r', required = True, metavar = 'string', help = 'The genome reference file used')
  parser.add_argument('--tools_directory', '-s', required = False, metavar = 'string', help = 'The path to the directory where the tools live')
  parser.add_argument('--mosaic_json', '-m', required = True, metavar = 'string', help = 'The json file describing the Mosaic parameters')

  # Optional mosaic arguments
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

# SpliceAI processing script version
version = "0.0.1"

if __name__ == "__main__":
  main()
