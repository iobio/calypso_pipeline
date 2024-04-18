from datetime import date
from os.path import exists
from sys import path

import argparse
import os
import sys
#import calypso_path as cpath
import calypso_resources as res

# Add the path of the common functions and import them
path.append(os.path.dirname(__file__) + '/components')
import tools_bcftools as bcftools
import calypso_mosaic_resources as mosr

def main():
  global allowedClasses
  global bcftoolsExe
  global impacts

  # Parse the command line
  args = parseCommandLine()

  # Define the executable bcftools command
  if args.tools_directory:
    if not args.tools_directory.endswith('/'): args.tools_directory = str(args.tools_directory) + '/'
    bcftoolsExe = args.tools_directory + 'bcftools/bcftools'
  else: bcftoolsExe = 'bcftools'

  # Open an output tsv file to write annotations to
  outputFile = open(args.output_tsv, 'w')

  # Loop over the annotations that are to be uploaded to Mosaic for this resource and get the annotation name, uid and
  # type
  annotations = {}
  annotations['consequence'] = {'uid': 'gene_consequence_GRCh38', 'type': 'string'}
  annotations['impact'] = {'uid': 'gene_impact_GRCh38', 'type': 'string'}

  # Write the header line to the tsv file
  print('CHROM\tSTART\tEND\tREF\tALT\t', '\t'.join(str(x['uid']) for x in annotations.values()), sep = '', file = outputFile)

  # Loop over all records in the vcf file
  for record in os.popen(bcftools.splitvepQueryGeneConsImp(bcftoolsExe, args.input_vcf)).readlines():

    # Split the record on tabs and remove the gene name
    fields = record.rstrip().split('\t')
    fields.pop(5)

    # If the consequence has multiple values, take the highest impact
    if '&' in fields[5]:
      consequences = fields[5].split('&')
      consequence  = consequences[0].lower()
      if consequence.endswith('_variant'):
        consequence = consequence[:-(len('_variant'))]
      outputConsequence = [consequence]
      if consequence not in impacts: fail('Unknown consequence: ' + str(consequence))
      impact = impacts[consequence]
      for option in consequences[1:]:
        option = option.lower()
        if option.endswith('_variant'):
          option = option[:-(len('_variant'))]
        if option not in impacts: fail('Unknown consequence: ' + str(option))
        if impacts[consequence] > impact:
          consequence = option
          impact      = impacts[option]
          outputConsequence = [option]
        elif impacts[consequence] == impact: outputConsequence.append(option)
      fields[5] = ','.join(outputConsequence)

    # If there is only a single value, strip '_variant'
    if fields[5].endswith('_variant'):
      fields[5] = fields[5][:-(len('_variant'))]

    # If all values are '.', this line can be ignored
    uniqueValues = set(fields[5:])
    if len(uniqueValues) == 1 and list(uniqueValues)[0] == '.': continue

    # Update the chromosome and position
    fields[0], fields[2] = updateCoords(fields[0], fields[2])
  
    # Build the output record from the updated fields
    print('\t'.join(fields), file = outputFile)

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
  parser.add_argument('--tools_directory', '-s', required = True, metavar = 'string', help = 'The path to the directory where the tools live')

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
  print(message, sep = '')
  exit(1)

# Initialise global variables

# Pipeline version
version = "0.0.1"

# Define the bcftools executable
bcftoolsExe = False

impacts = {}
impacts['transcript_ablation'] = 4
impacts['splice_acceptor'] = 4
impacts['splice_donor'] = 4
impacts['stop_gained'] = 4
impacts['frameshift'] = 4
impacts['stop_lost'] = 4
impacts['start_lost'] = 4
impacts['transcript_amplification'] = 4
impacts['feature_elongation'] = 4
impacts['feature_truncation'] = 4
impacts['inframe_insertion'] = 3
impacts['inframe_deletion'] = 3
impacts['missense'] = 3
impacts['protein_altering'] = 3
impacts['splice_donor_5th_base'] = 2
impacts['splice_region'] = 2
impacts['splice_donor_region'] = 2
impacts['splice_polypyrimidine_tract'] = 2
impacts['incomplete_terminal_codon'] = 2
impacts['start_retained'] = 2
impacts['stop_retained'] = 2
impacts['synonymous'] = 2
impacts['coding_sequence'] = 1
impacts['mature_mirna'] = 1
impacts['5_prime_utr'] = 1
impacts['3_prime_utr'] = 1
impacts['non_coding_transcript_exon'] = 1
impacts['intron'] = 1
impacts['nmd_transcript'] = 1
impacts['non_coding_transcript'] = 1
impacts['coding_transcript'] = 1
impacts['upstream_gene'] = 1
impacts['downstream_gene'] = 1
impacts['tfbs_ablation'] = 1
impacts['tfbs_amplification'] = 1
impacts['tf_binding_site'] = 1
impacts['regulatory_region_ablation'] = 1
impacts['regulatory_region_amplification'] = 1
impacts['regulatory_region'] = 1
impacts['intergenic'] = 1
impacts['sequence'] = 1

if __name__ == "__main__":
  main()
