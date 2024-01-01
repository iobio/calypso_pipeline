from os.path import exists

import argparse
import os
import json
import subprocess
import sys

# Add the path of the common functions and import them
from sys import path
path.append(os.path.dirname(__file__) + '/../components')
import tools_bcftools as bcftools
import calypso_path as cpath

def main():
  global mosaicConfig

  # Parse the command line
  args = parseCommandLine()

  # Build the path to include all components
  sys.path = cpath.buildPath(sys.path, args.utils_directory)
  import mosaic_config
  import api_genes as api_g
  import api_variant_annotations as api_va

  # Parse the Mosaic config file to get the token and url for the api calls
  mosaicRequired = {'MOSAIC_TOKEN': {'value': args.token, 'desc': 'An access token', 'long': '--token', 'short': '-t'},
                    'MOSAIC_URL': {'value': args.url, 'desc': 'The api url', 'long': '--url', 'short': '-u'}}
  mosaicConfig   = mosaic_config.mosaicConfigFile(args.config)
  mosaicConfig   = mosaic_config.commandLineArguments(mosaicConfig, mosaicRequired)

  # Define the executable bcftools command
  if not args.tools_directory.endswith('/'): args.tools_directory = str(args.tools_directory) + '/'
  bcftoolsExe = args.tools_directory + 'bcftools/bcftools'

  # Check thet the compound hets vcf file exists
  if not exists(args.input): fail('The comp hets vcf file ' + str(args.input) + ' was not found')

  # Check if a private annotation already exists. If so, delete it
  data = api_va.getAnnotations(mosaicConfig, args.project_id)
  for annotation in data:
    if annotation['privacy_level'] == 'private' and annotation['name'] == 'Compound Hets': api_va.deleteAnnotationById(mosaicConfig, args.project_id, annotation['id'])
  
  # Create a private annotation for storing comp het variants
  data = api_va.createPrivateAnnotationIdUid(mosaicConfig, 'Compound Hets', 'float', args.project_id)
  annotationUid = data['uid']
  annotationId  = data['id']

  # Open an output tsv file to include compound het annotation
  comphetsFilename = 'comphets.tsv'
  comphets         = open(comphetsFilename, 'w')

  # Write the header line to the tsv file
  print('CHROM\tSTART\tEND\tREF\tALT\t', annotationUid, sep = '', file = comphets)

  # Loop over all the variants in the file
  genes = {}
  for line in os.popen(bcftools.query(bcftoolsExe, args.input, ['BCSQ'])):
    fields = line.rstrip().split('\t')

    # Update the chromosome and position
    fields[0], fields[2] = updateCoords(fields[0], fields[2])

    # Get the gene to add to the gene set
    gene = fields[5].split('|')[1]
    if gene not in genes: genes[gene] = 1
    else: genes[gene] += 1

    # Update the INFO field to 1. All variants that are part of a comp het will have a value of 1, all
    # other variants will have no value
    fields[5] = '1'

    # Print the fields to the output tsv file
    print('\t'.join(fields), file = comphets)

  # Loop over the genes and ensure that each gene appears as least twice as it must to be a comp het
  failedGenes = []
  for gene in genes:
    if genes[gene] < 2: failedGenes.append(gene)
  if failedGenes: fail('The following genes only included a single variant so cannot be a comp het. Please check the comp hets vcf\n  ' + ', '.join(failedGenes))

  # Check if a gene set already exists. If so, delete it
  for geneSet in api_g.getGeneSets(mosaicConfig, args.project_id):
    if geneSet['name'] == 'Compound het genes': api_g.deleteGeneSet(mosaicConfig, args.project_id, geneSet['id'])

  # Create a gene set of the comp het genes
  description = 'Genes containing compound heterozygotes'
  geneSetId, skippedGenes = api_g.postGeneSetByName(mosaicConfig, args.project_id, 'Compound het genes', description, genes.keys(), 'Homo sapiens')

  # Close the annotations file
  comphets.close()

  # Upload the annotations to Mosaic
  api_va.uploadAnnotations(mosaicConfig, args.project_id, os.getcwd() + '/' + str(comphetsFilename), 'true')

# Input options
def parseCommandLine():
  global version

  parser = argparse.ArgumentParser(description='Process the command line')

  # Mosaic arguments
  parser.add_argument('--config', '-c', required = False, metavar = "string", help = "The config file for Mosaic")
  parser.add_argument('--token', '-t', required = False, metavar = "string", help = "The Mosaic authorization token")
  parser.add_argument('--url', '-u', required = False, metavar = "string", help = "The base url for Mosaic")

  # Required arguments
  parser.add_argument('--input', '-i', required = True, metavar = 'string', help = 'The exomiser variants.tsv file')
  parser.add_argument('--project_id', '-p', required = True, metavar = 'float', help = 'The Mosaic project id that these exomiser results will be uploaded to')

  # The public-utils and tools directories
  parser.add_argument('--utils_directory', '-l', required = True, metavar = 'string', help = 'The path to the public-utils directory')
  parser.add_argument('--tools_directory', '-s', required = True, metavar = 'string', help = 'The path to the tools directory')

  # Version
  parser.add_argument('--version', '-v', action='version', version='Exomiser variant output parser version: ' + str(version))

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

# Version
version = 0.01

# Store the ids and uids of the exomiser annotations
annIds = {}

# Store information related to Mosaic
mosaicConfig = {}

if __name__ == "__main__":
  main()
