#!/usr/bin/python

from __future__ import print_function
from datetime import date
from os.path import exists
from sys import path

import argparse
import os
import sys

# Add the path of the common functions and import them
path.append(os.path.dirname(__file__) + '/components')
import tools_bcftools as bcftools
import calypso_mosaic_resources as mosr

def main():
  global allowedClasses

  # Parse the command line
  args = parseCommandLine()

  # Read the mosaicJson file to get information on how to process different annotations
  mosaicInfo = mosr.readMosaicJson(args.mosaic_json, args.reference)

  # Check that the provided annotation class is valid
  processClass = mosaicInfo['resources'][args.resource]['class']
  if processClass in allowedClasses:

    # Open an output tsv file to write annotations to
    outputFile = open(args.output_tsv, 'w')
  
    # Write the header line to the tsv file
    print('CHROM\tSTART\tEND\tREF\tALT\t', '\t'.join(args.uids.replace(' ', '').split(',')), sep = '', file = outputFile)
  
    # Loop over the vcf and process according to the annotation class
    if processClass == "A": processClassA(args.input_vcf, args.tags.replace(' ', '').split(','), outputFile)
    elif processClass == "B": processClassB(args.input_vcf, args.tags.replace(' ', '').split(','), outputFile)
    elif processClass == "OMIM": processClassOMIM(args.input_vcf, args.tags.replace(' ', '').split(','), outputFile)
    elif processClass == "hgvs": processClassHgvs(mosaicInfo['resources'][args.resource], args.input_vcf, args.tags.replace(' ', '').split(','), outputFile)
  
    # Close the output tsv file
    outputFile.close()

  # Write a warning if the annotation class is not recognised and do not process the vcf file
  else: print('Unable to process annotations for resource: ' + args.resource, sep = '', file = sys.stderr)

# Input options
def parseCommandLine():
  global version
  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  parser.add_argument('--config', '-c', required = True, metavar = "string", help = "The config file for Mosaic")
  parser.add_argument('--input_vcf', '-i', required = True, metavar = 'string', help = 'The input vcf file to annotate')
  parser.add_argument('--output_tsv', '-o', required = True, metavar = 'string', help = 'The output tsv file')
  parser.add_argument('--resource', '-e', required = True, metavar = 'string', help = 'The name of the resource (used for the output tsv name')
  parser.add_argument('--reference', '-r', required = True, metavar = 'string', help = 'The genome reference file used')
  parser.add_argument('--tags', '-g', required = True, metavar = 'string', help = 'A comma separated list of VCF INFO tags to be extracted')
  parser.add_argument('--uids', '-d', required = True, metavar = 'string', help = 'A comma separated list of uids for the resource annotations')
  parser.add_argument('--mosaic_json', '-m', required = True, metavar = 'string', help = 'The json file describing the Mosaic parameters')

  # Optional mosaic arguments
  parser.add_argument('--token', '-t', required = False, metavar = "string", help = "The Mosaic authorization token")
  parser.add_argument('--url', '-u', required = False, metavar = "string", help = "The base url for Mosaic")
  parser.add_argument('--attributes_project', '-a', required = False, metavar = "integer", help = "The Mosaic project id that contains public attributes")

  # Version
  parser.add_argument('--version', '-v', action="version", version='Calypso annotation pipeline version: ' + str(version))

  return parser.parse_args()

# Process class A annotations. This is for annotations that are floats and if multiple values occur, only the maximum value should be
# uploaded to Mosaic. This includes the following annoatations:
#   - CCR
#   - EVE
#   - gnomAD
#   - MutScore
#   - pLI
#   - REVEL
def processClassA(vcf, tags, outputFile):

  # Loop over all records in the vcf file
  for record in os.popen(bcftools.query(vcf, tags)).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')
  
    # Check that the annotation has a value. If not, skip this record. The only records that will be output are those with values
    hasValue = False
    for i in range(5, len(fields)):
      if fields[5] != '.':
        hasValue = True
        break
    if hasValue:
  
      # Update the chromosome and position
      fields[0], fields[2] = updateCoords(fields[0], fields[2])
  
      # If there are multiple values, only output the largest
      for i in range(5, len(fields)):
        if ',' in fields[i]:
          outputValue = 0.
          for value in fields[i].split(','):
            if float(value) > float(outputValue): outputValue = value
          fields[i] = str(outputValue)
  
      # Check that the value is a float
      try: typeTest = float(fields[5])
      except: fail('Invalid valid for annotation: ' + record.rstrip())
  
      # Make sure the value falls between 1E-37 and 1E+37
      if float(fields[5]) <= 1e-37: fields[5] = "1e-37"
      elif float(fields[5]) >= 1e37: fields[5] = "1e37"
  
      # Build the output record from the updated fields
      print('\t'.join(fields), file = outputFile)

# Process class B annotations. This is for annotations that are strings that do not undergo any modifications
#   - dbSNP
def processClassB(vcf, tags, outputFile):

  # Loop over all records in the vcf file
  for record in os.popen(bcftools.query(vcf, tags)).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')
  
    # Check that the variant has a value. If not, skip this record. The only records that will be output are those with values
    hasValue = False
    for i in range(5, len(fields)):
      if fields[5] != '.':
        hasValue = True
        break
    if hasValue:
  
      # Update the chromosome and position
      fields[0], fields[2] = updateCoords(fields[0], fields[2])
  
      # Build the output record from the updated fields
      print('\t'.join(fields), file = outputFile)

# Process OMIM annotations
def processClassOMIM(vcf, tags, outputFile):

  # Loop over all records in the vcf file
  for record in os.popen(bcftools.query(vcf, tags)).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')
  
    # Check that the variant has a value. If not, skip this record. The only records that will be output are those with values
    hasValue = False
    for i in range(5, len(fields)):
      if fields[5] != '.':
        hasValue = True
        break
    if hasValue:
  
      # Update the chromosome and position
      fields[0], fields[2] = updateCoords(fields[0], fields[2])
  
      # If the field begins with a comma, or has two consecutive commas, there is no value. For example, a variant may be associated
      # with three MIM ids, but only the third has an associated inheritance. The inheritance will appear as ',,AR', for example. In
      # this case, the leading ',' should be replaced with a leading 'None', and the consecutive commas should be replaced with
      # ',None,' yielding a value of None,None,AR for upload to Mosaic
      for i in range(5, len(fields)):
        if fields[i].startswith(','): fields[i] = 'None' + fields[i]
        fields[i] = fields[i].replace(',,', ',None,')
  
      # Build the output record from the updated fields
      print('\t'.join(fields), file = outputFile)

# Process HGVS annotations
def processClassHgvs(resourceInfo, vcf, tags, outputFile):

  # Loop over all records in the vcf file
  for record in os.popen(bcftools.splitvepQuery(vcf, [resourceInfo['info_field']])).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')
  
    # Check that the variant has a value. If not, skip this record. The only records that will be output are those with values
    if fields[5] != '.':
  
      # Update the chromosome and position
      fields[0], fields[2] = updateCoords(fields[0], fields[2])
      annotations = fields.pop().split('|')
      hasValue = False
      for tag in tags:
        value = annotations[resourceInfo['annotations'][tag]['position'] - 1]
        fields.append(annotations[resourceInfo['annotations'][tag]['position'] - 1])
        if value: hasValue = True
  
      # Build the output record from the updated fields
      if hasValue: print('\t'.join(fields), file = outputFile)

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

# Pipeline version
version = "1.0.0"

# Define the allowed annotation classes
allowedClasses = ["A", "B", "OMIM", "hgvs"]

if __name__ == "__main__":
  main()
