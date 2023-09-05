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
path.append(os.path.dirname(__file__) + '/../components')
import calypso_path as cpath
import tools_bcftools as bcftools

def main():
  global mosaicConfig
  global patientTerms

  # Parse the command line
  args = parseCommandLine()

  # Build the path to include all components
  sys.path = cpath.buildPath(sys.path, args.utils_directory)
  import mosaic_config
  import api_sample_hpo_terms as api_hpo
  import api_samples as api_s
  import api_sample_attributes as api_sa
  import api_variant_annotations as api_va

  # Define the executable bcftools command
  if args.tools_directory:
    if not args.tools_directory.endswith('/'): args.tools_directory = str(args.tools_directory) + '/'
    bcftoolsExe = args.tools_directory + 'bcftools/bcftools'
  else: bcftoolsExe = 'bcftools'

  # Parse the Mosaic config file to get the token and url for the api calls
  mosaicRequired = {'MOSAIC_TOKEN': {'value': args.token, 'desc': 'An access token', 'long': '--token', 'short': '-t'},
                    'MOSAIC_URL': {'value': args.url, 'desc': 'The api url', 'long': '--url', 'short': '-u'}}
  mosaicConfig   = mosaic_config.mosaicConfigFile(args.config)
  mosaicConfig   = mosaic_config.commandLineArguments(mosaicConfig, mosaicRequired)

  # Get all the annotations in the project and find the HPO annotations. Store the uids of these
  annotations = api_va.getAnnotationDictNameIdUid(mosaicConfig, args.project_id)
  labelsUid   = False
  overlapsUid = False
  termsUid    = False
  for annotation in annotations:
    if annotation == 'HPO Labels': labelsUid = annotations[annotation]['uid']
    if annotation == 'HPO Terms': termsUid = annotations[annotation]['uid']
    if annotation == 'HPO Overlaps': overlapsUid = annotations[annotation]['uid']
  if not labelsUid: fail('No private annotation with the name "HPO Labels" present in this project')
  if not overlapsUid: fail('No private annotation with the name "HPO Overlaps" present in this project')
  if not termsUid: fail('No private annotation with the name "HPO Terms" present in this project')

  # Get the HPO terms for the proband
  sampleAttributes = api_sa.getSampleAttributesDictNameId(mosaicConfig, args.project_id)
  if 'Relation' not in sampleAttributes: fail('The sample attribute "Relation" is not present in this project. This attribute must be present and populated')
  relationId = sampleAttributes['Relation']

  # Find the proband as the sample with the 
  sampleIds = api_s.getSamplesDictIdName(mosaicConfig, args.project_id)
  probandId = False
  for sampleId in sampleIds:
    sampleAttributes = api_sa.getAttributesForSampleDictId(mosaicConfig, args.project_id, sampleId)
    if relationId not in sampleAttributes: fail('The sample attribute "Relation" is not present for sample ' + str(sampleIds[sampleId]))
    else:
      if sampleAttributes[relationId]['values'][0]['value'] == 'Proband':
        probandId = sampleId
        break
  if not probandId: fail('Could not determine the proband from the "Relation" sample attribute')

  # Get the HPO terms for the proband
  patientTerms = api_hpo.getSampleHpoList(mosaicConfig, args.project_id, probandId)
  if not patientTerms: fail('No HPO terms are associated with the sample')

  # Get HPO / gene associations
  parseHpoGene(args.hpo)

  # Parse vcf file
  parseVcf(args.input_vcf, args.output_tsv, bcftoolsExe, labelsUid, overlapsUid, termsUid)

# Input options
def parseCommandLine():
  global version

  parser = argparse.ArgumentParser(description='Process the command line')

  # Mosaic arguments
  parser.add_argument('--config', '-c', required = False, metavar = "string", help = "The config file for Mosaic")
  parser.add_argument('--token', '-t', required = False, metavar = "string", help = "The Mosaic authorization token")
  parser.add_argument('--url', '-u', required = False, metavar = "string", help = "The base url for Mosaic")

  # The project id
  parser.add_argument('--project_id', '-p', required = True, metavar = "string", help = "The project id that variants will be uploaded to")

  # Directories where required tools live
  parser.add_argument('--utils_directory', '-l', required = True, metavar = 'string', help = 'The path to the public-utils directory')
  parser.add_argument('--tools_directory', '-s', required = False, metavar = 'string', help = 'The path to the tools directory')

  # The file containing HPO:gene associations
  parser.add_argument('--hpo', '-f', required = True, metavar = "string", help = "The gene:HPO associations file")

  # The input vcf to base annotations on
  parser.add_argument('--input_vcf', '-i', required = True, metavar = 'string', help = 'The filtered vcf file')

  # The output tsv for upload to Mosaic
  parser.add_argument('--output_tsv', '-o', required = True, metavar = 'string', help = 'The output tsv file')

  # Version
  parser.add_argument('--version', '-v', action="version", version='HPO prioritization pipeline version: ' + str(version))

  return parser.parse_args()

# Get HPO / gene associations
def parseHpoGene(hpo):
  global hpoGene
  global hpoLabels
  global patientTerms

  # Open the file containing HPO:gene associations
  hpoFile = open(hpo, 'r')
  while True:
    line = hpoFile.readline().rstrip()
    if not line: break

    # Ignore the header
    if line.startswith('#'): continue
    hpo = line.split('\t')[2]

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
def parseVcf(vcf, tsv, bcftoolsExe, labelsUid, overlapsUid, termsUid):
  global hpoGene
  global hpoLabels
  global patientTerms
  global samples

  # Open output files. There needs to be one for public and one for private annotations
  hpoFile  = open(tsv, 'w')

  # Write the header lines to the output files
  headerBase = 'CHROM\tSTART\tEND\tREF\tALT'
  print(headerBase, termsUid, labelsUid, overlapsUid, sep = '\t', file = hpoFile)

  # Read the standard input
  for line in os.popen(bcftools.query(bcftoolsExe, vcf, ['BCSQ'])).readlines():

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
