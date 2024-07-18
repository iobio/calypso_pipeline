from datetime import date
from os.path import exists
from sys import path
from pprint import pprint

import argparse
import os
import sys

import read_resource_jsons as read_resources

def main():

  # Parse the command line
  args = parseCommandLine()

  # Define the executable bcftools command
  if args.tools_directory:
    if not args.tools_directory.endswith('/'): args.tools_directory = str(args.tools_directory) + '/'
    bcftools = args.tools_directory + 'bcftools/bcftools'
  else:
    bcftools = 'bcftools'

  # Open an output tsv file to write annotations to
  outputFile = open('variant_quality.tsv', 'w')

  # Write the header line to the tsv file
  print('CHROM\tSTART\tEND\tREF\tALT\t', args.uid, file = outputFile)

  # Loop over all records in the vcf file

  ### UPDATE FOR ONLY HIGH AND LOW QUAL - NO MEDIUM
  #command = bcftools + ' query -f \'%CHROM\\t%POS\\t%END\\t%REF\\t%ALT\\t%INFO/het_low_qual\\t%INFO/het_med_qual\\t%INFO/het_hi_qual\\t%INFO/hom_low_qual\\t%INFO/hom_med_qual\\t%INFO/hom_hi_qual' + '\\n\' ' + str(args.input_vcf)
  command = bcftools + ' query -f \'%CHROM\\t%POS\\t%END\\t%REF\\t%ALT\\t%INFO/het_low_qual\\t%INFO/het_hi_qual\\t%INFO/hom_low_qual\\t%INFO/hom_hi_qual' + '\\n\' ' + str(args.input_vcf)
  for record in os.popen(command).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')
    fields[0], fields[2] = updateCoords(fields[0], fields[2])

    # Since all samples in the vcf are evaluated for quailty, there could be sample ids associated with multiple of
    # the variant quality annotations
    variant_quality = 'Fail'
    if args.proband in fields[5] or args.proband in fields[7]:
      variant_quality = 'Fail'
    elif args.proband in fields[6] or args.proband in fields[8]:
      variant_quality = 'Pass'

    # If the proband quality is fail, but there are other samples that are high quality, mark this as Fail+
    if variant_quality == 'Fail':
      if fields[6] != '.' or fields[8] != '.':
        variant_quality = 'Fail+'
    ### UPDATE FOR ONLY HIGH AND LOW QUAL - NO MEDIUM
    #text = ''.join(fields[5:11])
    #if text.count('.') != 5:
    #if text.count('.') != 3:
    #  fail('Unexpected number of variant confidence values assigned to variant at ' + str(fields[0]) + ':' + str(fields[1]))

    # If there are no .'s at the beginning of "text", this is a het_low_qual. If there is one, then this is a het_med_qual etc
    # based on the order in the command above this loop. All samples in the vcf are evaluated for quality, so there could be a
    # comma separated list of samples associated with the quality tag
    #variant_quality = False

    # THE FOLLOWING IS FOR THE CASE WHERE WE HAVE MEDIUM
#    if text.startswith('.....'):
#      variant_quality = 'High'
#    elif text.startswith('....'):
#      variant_quality = 'Medium'
#    elif text.startswith('...'):
#      variant_quality = 'Low'
#    elif text.startswith('..'):
#      variant_quality = 'High'
#    elif text.startswith('.'):
#      variant_quality = 'Medium'
#    else:
#      variant_quality = 'Low'

    # AND THE FOLLOWING IS JUST PASS \ FAIL
#    if text.startswith('...'):
#      variant_quality = 'Pass'
#    elif text.startswith('..'):
#      variant_quality = 'Fail'
#    elif text.startswith('.'):
#      variant_quality = 'Pass'
#    else:
#      variant_quality = 'Fail'
  
    print('\t'.join(fields[0:5]), '\t', variant_quality, sep = '', file = outputFile)

  # Close the output tsv file
  outputFile.close()

# Input options
def parseCommandLine():
  global version
  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  parser.add_argument('--input_vcf', '-i', required = True, metavar = 'string', help = 'The input vcf file to annotate')
  parser.add_argument('--tools_directory', '-t', required = True, metavar = 'string', help = 'The path to the directory where the tools live')
  parser.add_argument('--uid', '-u', required = True, metavar = 'string', help = 'The uid of the variant quality attribute')
  parser.add_argument('--proband', '-p', required = True, metavar = 'string', help = 'The name of the proband as it appears in the vcf header')

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

if __name__ == "__main__":
  main()
