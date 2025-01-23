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

  # Open output tsv files to write annotations to
  vq_output_file = open('variant_quality.tsv', 'w')
  print('CHROM\tSTART\tEND\tREF\tALT\t', args.uid, file = vq_output_file)

  # Store the samples that must have the minimum GQ value
  if args.parents:
    parents = args.parents.split(',') if ',' in args.parents else [args.parents]
    trio_output_file = open('trio_quality.tsv', 'w')
    if not args.trio_uid:
      fail('parents have been provided and so --trio_uid (-d) must also be provided')
    print('CHROM\tSTART\tEND\tREF\tALT\t', args.trio_uid, file = trio_output_file)
  if args.family:
    family = args.family.split(',') if ',' in args.family else [args.parents]
    family_output_file = open('family_quality.tsv', 'w')
    if not args.family_uid:
      fail('additional family members have been provided and so the --family_uid (-m) must also be provided')
    print('CHROM\tSTART\tEND\tREF\tALT\t', args.family_uid, file = family_output_file)

  # Loop over all records in the vcf file
  command = bcftools + ' query -f \'%CHROM\\t%POS\\t%END\\t%REF\\t%ALT\\t%INFO/het_low_qual\\t%INFO/het_hi_qual\\t%INFO/hom_low_qual\\t%INFO/hom_hi_qual\\t%INFO/fam_geno' + '\\n\' ' + str(args.input_vcf)
  for record in os.popen(command).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')
    fields[0], fields[2] = updateCoords(fields[0], fields[2])

    # Since all samples in the vcf are evaluated for quailty, there could be sample ids associated with multiple of
    # the variant quality annotations
    variant_quality = 'Fail'
    trio_quality = False
    family_quality = False
    if args.proband in fields[5] or args.proband in fields[7]:
      variant_quality = 'Fail'
    elif args.proband in fields[6] or args.proband in fields[8]:
      variant_quality = 'Pass'

      # If additional samples were provided, check if they are present in the fam_geno annotation. If the proband quality is
      # Fail, there is no need for this check
      if args.parents:
        trio_quality = True
        for sample in parents:
          if sample not in fields[9]:
            trio_quality = False
            break
      if args.family:
        family_quality = True
        for sample in family:
          if sample not in fields[9]:
            family_quality = False
            break

    # If the proband quality is fail, but there are other samples that are high quality, mark this as Fail+
    if variant_quality == 'Fail':
      if fields[6] != '.' or fields[8] != '.':
        variant_quality = 'Fail+'
  
    print('\t'.join(fields[0:5]), '\t', variant_quality, sep = '', file = vq_output_file)
    if trio_quality:
      print('\t'.join(fields[0:5]), '\tPass', sep = '', file = trio_output_file)
    if family_quality:
      print('\t'.join(fields[0:5]), '\tPass', sep = '', file = family_output_file)

  # Close the output tsv file
  vq_output_file.close()
  if args.parents:
    trio_output_file.close()
  if args.family:
    family_output_file.close()

# Input options
def parseCommandLine():
  global version
  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  parser.add_argument('--input_vcf', '-i', required = True, metavar = 'string', help = 'The input vcf file to annotate')
  parser.add_argument('--tools_directory', '-t', required = True, metavar = 'string', help = 'The path to the directory where the tools live')
  parser.add_argument('--proband', '-p', required = True, metavar = 'string', help = 'The name of the proband as it appears in the vcf header')

  # The uids of the quality annotations. The trio and family quality uids should only be included if --parents and --family are used
  parser.add_argument('--uid', '-u', required = True, metavar = 'string', help = 'The uid of the variant quality attribute')
  parser.add_argument('--trio_uid', '-d', required = False, metavar = 'string', help = 'The uid of the trio variant quality attribute')
  parser.add_argument('--family_uid', '-m', required = False, metavar = 'string', help = 'The uid of the family variant quality attribute')

  # Optionally list the sample names of the mother and father. If these are provided, generate the additional family quality annotation to
  # determine whether the inheritance can be believed
  parser.add_argument('--parents', '-a', required = False, metavar = 'string', help = 'A comma separated list of the sample names of the mother and father')
  parser.add_argument('--family', '-f', required = False, metavar = 'string', help = 'A comma separated list of the sample names of additional family members that should be considered in the family quality filter')

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
  print('ERROR: ', message, sep = '')
  exit(1)

# Initialise global variables

if __name__ == "__main__":
  main()
