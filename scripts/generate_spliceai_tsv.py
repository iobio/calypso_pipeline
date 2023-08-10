#!/usr/bin/python

from __future__ import print_function
from datetime import date
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
  annotations         = []
  i = 1
  for annotation in mosaicInfo['resources']['SpliceAI']['annotations']:
    info = {'name': annotation, 'position': False, 'uid': mosaicInfo['resources']['SpliceAI']['annotations'][annotation]['uid'], 'operation': False}

    # If this annotation is to be constructed, determine the instructions and add them to this annotation
    operation = mosaicInfo['resources']['SpliceAI']['annotations'][annotation]['operation']
    if operation: info['operation'] = operation

    # Otherwise, this is an annotation to be included as is
    else:
      annotationList.append(annotation)
      annotationPositions[annotation] = int(i) + 4
      info['position'] = int(i) + 4
      i += 1

    # Update the string of annotation uids and store the info on this annotation
    annotationUids += str(mosaicInfo['resources']['SpliceAI']['annotations'][annotation]['uid']) + '\t'
    annotations.append(info)

  # Strip the trailing comma from the string listing the uids of all annotations
  annotationUids = annotationUids.rstrip('\t')

  # Loop over the annotations and process those that require an "operation", e.g. the max of other annotations
  for annotation in annotations:
    if annotation['operation']:

      # For the "max" operation, there should be a list (in "fields"), that indicates which annotations to find the max of.
      # Find the positions of these fields in annotationList and these will be used to extract the values from the bcftools
      # query record
      if annotation['operation'] == 'max': maxInstructions(annotation, mosaicInfo['resources']['SpliceAI']['annotations'][annotation['name']], annotationPositions)

      # Fail if the operation is not known
      else: fail('Unknown operation (' + operation + ') for annotation ' + annotation)

  # Open an output tsv file to write annotations to
  outputFile = open(args.output_tsv, 'w')
  
  # Write the header line to the tsv file
  print('CHROM\tSTART\tEND\tREF\tALT\t', annotationUids, sep = '', file = outputFile)

  # Loop over all records in the vcf file
  for record in os.popen(bcftools.query(bcftoolsExe, args.input_vcf, annotationList)).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')

    # Update the chromosome and position. If the chromosome begins with "chr", this needs to be stripped off. The end coordinate
    # (fields[2] - the start coordinate is fields[1]) needs to have one added to create a 0-based, half open coordinate
    fields[0], fields[2] = updateCoords(fields[0], fields[2])

    # The output line for each vcf record is the first 5 fields (coordinates and alleles), then the annotations in order.
    # Create the start of this line
    output = '\t'.join(fields[0:5])

    # If all values are '.', this line can be ignored
    uniqueValues = set(fields[5:])
    if len(uniqueValues) == 1 and list(uniqueValues)[0] == '.': continue

    # Loop over the annotations in the order they need to be output. For annotations that just need the value posting,
    # "position" will be set and defines the position in the fields array. If not set, this is an annotation that needs
    # to be generated, so this will be performed
    for annotation in annotations:
      if annotation['position']:

        # If multiple values are available for this position, the values for this position will be separated by a comma
        # In this scenario, take the maximum value
        if ',' in fields[annotation['position']]:
          maxValue = 0
          for value in fields[annotation['position']].split(','):
            if abs(float(value)) > abs(maxValue): maxValue = float(value)
          output += '\t' + str(maxValue)
        else: output += '\t' + str(fields[annotation['position']])

      # If multiple values are available for this position, the values for this position will be separated by a comma
      # In this scenario, take the maximum value of the maximum values
      elif annotation['operation'] == 'max':
        maxValue = 0
        for position in annotation['fields']:
          if ',' in fields[position]:
            for value in fields[position].split(','):
              if abs(float(value)) > maxValue: maxValue = float(value)
          else:
            if abs(float(fields[position])) > maxValue: maxValue = float(fields[position])
        output += '\t' + str(maxValue)
      else: fail('Annotation ' + annotation['name'] + 'does not have a position set, nor has a known operation to perform')

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

def maxInstructions(annotation, info, positions):
  fields = []
  if 'fields' not in info: fail('No "fields" section for annotation ' + str(annotation) + '. These are required to perform the max operation')
  for field in info['fields']: fields.append(positions[field])
  annotation['fields'] = fields

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
