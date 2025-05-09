from datetime import date
from os.path import exists
from sys import path

import argparse
import os
import sys

import read_resource_jsons as read_resources

def main():

  # Parse the command line
  args = parseCommandLine()

  # Define the executable bcftools command
  if args.tools_directory:
    if not args.tools_directory.endswith('/'):
      args.tools_directory = str(args.tools_directory) + '/'
    bcftools = args.tools_directory + 'bcftools/bcftools'
  else:
    bcftools = 'bcftools'

  # Read the mosaicJson file to get information on how to process different annotations
  mosaicInfo = read_resources.read_mosaic_json(args.mosaic_json, args.reference)

  # Open an output tsv file to write annotations to
  output_file = open(args.output_tsv, 'w')

  # Loop over the annotations that are to be uploaded to Mosaic for this resource and get the annotation name, uid and
  # type
  annotations = {}
  uids        = []
  for annotation in mosaicInfo['resources'][args.resource]['annotations']:
    uid     = mosaicInfo['resources'][args.resource]['annotations'][annotation]['uid']
    tagType = mosaicInfo['resources'][args.resource]['annotations'][annotation]['type']
    annotations[annotation] = {'uid': uid, 'type': tagType}

  # Write the header line to the tsv file
  print('CHROM\tSTART\tEND\tREF\tALT\t', '\t'.join(str(x['uid']) for x in annotations.values()), sep = '', file = output_file)

  # Parse the vcf file, check the annotations and generate a tsv file for upload to Mosaic
  process_vcf(bcftools, args.input_vcf, annotations, output_file)

  # Close the output tsv file
  output_file.close()

# Input options
def parseCommandLine():
  global version
  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  parser.add_argument('--input_vcf', '-i', required = True, metavar = 'string', help = 'The input vcf file to annotate')
  parser.add_argument('--output_tsv', '-o', required = True, metavar = 'string', help = 'The output tsv file')
  parser.add_argument('--resource', '-e', required = True, metavar = 'string', help = 'The name of the resource (used for the output tsv name')
  parser.add_argument('--reference', '-r', required = True, metavar = 'string', help = 'The genome reference file used')
  parser.add_argument('--tools_directory', '-s', required = True, metavar = 'string', help = 'The path to the directory where the tools live')
  parser.add_argument('--mosaic_json', '-m', required = True, metavar = 'string', help = 'The json file describing the Mosaic parameters')

  return parser.parse_args()

# The majority of annotations can be processed in a standard way. The type is checked and numerical
# annotations are checked to ensure they fall within the required bounds
def process_vcf(bcftools, vcf, annotations, output_file):

  # Loop over all records in the vcf file
  command = bcftools + ' query -f \'%CHROM\\t%POS\\t%END\\t%REF\\t%ALT\\t%INFO/' + '\\t%INFO/'.join(annotations.keys()) + '\\n\' ' + str(vcf)
  #print(command)
  for record in os.popen(command).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')

    # If all values are '.', this line can be ignored
    uniqueValues = set(fields[5:])
    if len(uniqueValues) == 1 and list(uniqueValues)[0] == '.':
      continue

    # Update the chromosome and position
    fields[0], fields[2] = updateCoords(fields[0], fields[2])
  
    # Loop over the annotations in the order they appear in the bcftools query command
    i = 5
    for annotation in annotations:

      # If the annotation is numeric
      if str(annotations[annotation]['type']) == 'integer' or str(annotations[annotation]['type']) == 'float':

        # If the are multiple values (e.g. there is a ',') in the value, ensure that all values are within
        # the limits (1E-37 < x < 1E37), but output all values
        values = fields[i].split(',') if ',' in fields[i] else [fields[i]]
        finalValues = []
        for value in values:
          if value == '.':
            finalValues.append('')
          else:
            try:
              typeTest = float(value)
            except:
              fail('Invalid value for annotation: ' + record.rstrip())
            if abs(float(value)) <= 1e-37:
              finalValues.append('0')
            elif abs(float(value)) >= 1e37:
              finalValues.append('1e37')
            else:
              finalValues.append(value)

        # Join the values back together with commas
        fields[i] = ','.join(finalValues)

      # Else if the annotation is a string, the annotation should be output as is
      elif str(annotations[annotation]['type']) == 'string':
        pass

      # If the type is unknown, fail
      else:
        fail('Annotation type (' + str(annotations[annotation]['type']) + ') for annotation ' + str(annotation) + ', is unknown')

      # Iterate to the next annotation field
      i += 1

    # Build the output record from the updated fields
    print('\t'.join(fields), file = output_file)

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

if __name__ == "__main__":
  main()
