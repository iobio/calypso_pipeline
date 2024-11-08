from datetime import date
from os.path import exists
from sys import path

import argparse
import os
import sys

import read_resource_jsons as read_resources

def main():

  # Parse the command line
  args = parse_command_line()

  # Define the executable bcftools command
  if args.tools_directory:
    if not args.tools_directory.endswith('/'):
      args.tools_directory = str(args.tools_directory) + '/'
    bcftools = args.tools_directory + 'bcftools/bcftools'
  else:
    bcftools = 'bcftools'

  # Open an output tsv file to write annotations to
  output_file = open(args.output_tsv, 'w')

  # Write the header line to the tsv file
  print('CHROM\tSTART\tEND\tREF\tALT\t', str(args.uid), sep = '', file = output_file)

  # Parse the vcf file, check the annotations and generate a tsv file for upload to Mosaic
  process_vcf(bcftools, args.input_vcf, output_file)

  # Close the output tsv file
  output_file.close()

# Input options
def parse_command_line():
  global version
  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  parser.add_argument('--input_vcf', '-i', required = True, metavar = 'string', help = 'The input vcf file to annotate')
  parser.add_argument('--output_tsv', '-o', required = True, metavar = 'string', help = 'The output tsv file')
  parser.add_argument('--tools_directory', '-s', required = True, metavar = 'string', help = 'The path to the directory where the tools live')
  parser.add_argument('--uid', '-u', required = True, metavar = 'string', help = 'The uid for the private comp het annotation')

  return parser.parse_args()

# The majority of annotations can be processed in a standard way. The type is checked and numerical
# annotations are checked to ensure they fall within the required bounds
def process_vcf(bcftools, vcf, output_file):

  # Loop over all records in the vcf file
  command = bcftools + ' query -f \'%CHROM\\t%POS\\t%END\\t%REF\\t%ALT\\tPass\\n\' ' + str(vcf)
  for record in os.popen(command).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')

    # Update the chromosome and position
    fields[0], fields[2] = updateCoords(fields[0], fields[2])
  
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
