#!/usr/bin/python

from __future__ import print_function
from datetime import date
from os.path import exists

import sys
import os
import argparse
import json
import math
import glob

# Add the path of the common functions and import them
from sys import path
path.append(os.path.dirname(__file__) + '/components')
import calypso_path as cpath
import tools_bcftools as bcftools

def main():
  global mosaicConfig
  global mosaicInfo

  # Parse the command line
  args = parseCommandLine()

  # Build the path to include all components
  sys.path = cpath.buildPath(sys.path, args.utils_directory)
  import mosaic_config
  import api_samples as api_s

  # Parse the Mosaic config file to get the token and url for the api calls
  mosaicRequired = {'MOSAIC_TOKEN': {'value': args.token, 'desc': 'An access token', 'long': '--token', 'short': '-t'},
                    'MOSAIC_URL': {'value': args.url, 'desc': 'The api url', 'long': '--url', 'short': '-u'}}
  mosaicConfig = mosaic_config.parseConfig(args.config, mosaicRequired)

  # Determine the family structure
  getProband(args, api_s)

  # Parse the patients HPO terms
  parsePatientHpoTerms(args.hpo_terms)

  # Get HPO / gene associations
  parseHpoGene(args)

  # Parse vcf file
  parseVcf(args)

# Input options
def parseCommandLine():
  global version

  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  parser.add_argument('--input_vcf', '-i', required = True, metavar = 'string', help = 'The filtered Calypso vcf file')
  parser.add_argument('--output_tsv', '-o', required = True, metavar = 'string', help = 'The output tsv file')
  parser.add_argument('--hpo', '-f', required = True, metavar = "string", help = "The gene:HPO associations file")
  parser.add_argument('--project_id', '-p', required = True, metavar = "string", help = "The project id that variants will be uploaded to")
  parser.add_argument('--hpo_terms', '-r', required = True, metavar = "string", help = "The HPO terms associated with the patient")
  parser.add_argument('--terms_uid', '-e', required = True, metavar = "string", help = "The uid for the HPO terms annotation")
  parser.add_argument('--labels_uid', '-b', required = True, metavar = "string", help = "The uid for the HPO labels annotation")
  parser.add_argument('--overlaps_uid', '-d', required = True, metavar = "string", help = "The uid for the HPO overlaps annotation")
  parser.add_argument('--utils_directory', '-l', required = True, metavar = 'string', help = 'The path to the public-utils directory')

  # Optional mosaic arguments
  parser.add_argument('--config', '-c', required = False, metavar = "string", help = "The config file for Mosaic")
  parser.add_argument('--token', '-t', required = False, metavar = "string", help = "The Mosaic authorization token")
  parser.add_argument('--url', '-u', required = False, metavar = "string", help = "The base url for Mosaic")

  # Version
  parser.add_argument('--version', '-v', action="version", version='HPO prioritization pipeline version: ' + str(version))

  return parser.parse_args()
  
# Parse the patients HPO terms
def parsePatientHpoTerms(hpoTerms):
  global patientTerms

  patientTerms = hpoTerms.split(',')

# Determine the id of the proband
def getProband(args, api_s):
  global mosaicConfig
  global samples

  # Get the samples Mosaic id and store
  samples = {}
  for sample in api_s.getSampleIds(mosaicConfig, args.project_id): samples[sample['name']] = sample['id']

# Get HPO / gene associations
def parseHpoGene(args):
  global hpoGene
  global hpoLabels
  global patientTerms

  # Open the file containing HPO:gene associations
  hpoFile = open(args.hpo, 'r')
  while True:
    line = hpoFile.readline().rstrip()
    if not line: break

    # Ignore the header
    if line.startswith('#'): continue
    hpo   = line.split('\t')[2]

    # Only store information on genes that are associated with HPO terms present in the patient terms
    if hpo in patientTerms:
      gene  = line.split('\t')[1]
      label = line.split('\t')[3]

      # Store the information
      if gene not in hpoGene: hpoGene[gene] = [hpo]
      else: hpoGene[gene].append(hpo)
      if hpo not in hpoLabels: hpoLabels[hpo] = label

  # Close the file
  hpoFile.close()

# Parse vcf file
def parseVcf(args):
  global hpoGene
  global hpoLabels
  global patientTerms
  global samples

  # Open output files. There needs to be one for public and one for private annotations
  hpoFile  = open(args.output_tsv, 'w')

  # Write the header lines to the output files
  headerBase = 'CHROM\tSTART\tEND\tREF\tALT'
  print(headerBase, args.terms_uid, args.labels_uid, args.overlaps_uid, sep = '\t', file = hpoFile)

  # Read the standard input
  for line in os.popen(bcftools.query(args.input_vcf, ['BCSQ'])).readlines():

    # Skip header lines
    if line.startswith('#'): continue
    line = line.rstrip()

    # Split the line up
    fields = line.split('\t')

    # If there is no BCSQ in the record, a '.' will be supplied and these should be ignored
    bcsq  = fields[5]
    genes = []
    if bcsq != '.':
      if ',' in bcsq:
        for anno in bcsq.split(','):

          # Get all the genes
          gene = anno.split('|')[1]
          if gene not in genes: genes.append(gene)
      else: genes = [bcsq.split('|')[1]]

      # Check if any of the genes have an association with an HPO term association with the patient
      hasAssociation    = False
      associatedTerms   = []
      associatedLabels  = []
      noAssociatedTerms = 0
      for gene in genes:
        if gene in hpoGene:
          for hpoTerm in hpoGene[gene]:
            if hpoTerm in patientTerms:
              hasAssociation = True
              if hpoTerm not in associatedTerms:
                noAssociatedTerms += 1
                associatedTerms.append(hpoTerm)
              if hpoLabels[hpoTerm] not in associatedLabels: associatedLabels.append(hpoLabels[hpoTerm])
  
      # Get the basic information required for upload to Mosaic as an annotation
      chrom = line.split('\t')[0]
      if chrom.startswith('chr'): chrom = chrom.strip('chr')
      start = line.split('\t')[1]
      end   = int(line.split('\t')[2]) + 1
      ref   = line.split('\t')[3]
      alt   = line.split('\t')[4]

      # If the label has more than 255 characters, trim it
      label = ','.join(associatedLabels)
      if len(label) > 255: label = label[0:252] + '...'
  
      # If there is an association with a patient HPO term, write out the tsv for upload to Mosaic
      if hasAssociation: print(chrom, start, end, ref, alt, ','.join(associatedTerms), label, noAssociatedTerms, sep = '\t', file = hpoFile)

  # Close output files
  hpoFile.close()

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = '')
  exit(1)

# Initialise global variables

# Version
version = "0.0.1"
date    = str(date.today())

# Store information related to Mosaic
mosaicConfig = {}
mosaicInfo   = {}

# Store the patients HPO terms
patientTerms = []

# Store the hpo:gene associations
hpoGene = {}

# Store the HPO labels
hpoLabels = {}

# Mosaic samples
samples = {}

if __name__ == "__main__":
  main()
