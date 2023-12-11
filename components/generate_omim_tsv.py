from datetime import date
from os.path import exists
from sys import path

import argparse
import os
import sys
import calypso_resources as res

# Add the path of the common functions and import them
path.append(os.path.dirname(__file__) + '/components')
import tools_bcftools as bcftools
import calypso_mosaic_resources as mosr
import calypso_path as cpath

def main():
  global allowedClasses
  global bcftoolsExe

  # Parse the command line
  args = parseCommandLine()

  # Define the executable bcftools command
  if args.tools_directory:
    if not args.tools_directory.endswith('/'): args.tools_directory = str(args.tools_directory) + '/'
    bcftoolsExe = args.tools_directory + 'bcftools/bcftools'
  else: bcftoolsExe = 'bcftools'
  sys.path = cpath.buildPath(sys.path, args.utils_directory)
  import mosaic_config
  import api_variant_annotations as api_va

  # Parse the Mosaic config file to get the token and url for the api calls
  mosaicRequired = {'MOSAIC_TOKEN': {'value': args.token, 'desc': 'An access token', 'long': '--token', 'short': '-t'},
                    'MOSAIC_URL': {'value': args.url, 'desc': 'The api url', 'long': '--url', 'short': '-u'},
                    'MOSAIC_ATTRIBUTES_PROJECT_ID': {'value': args.attributes_project, 'desc': 'The public attribtes project id', 'long': '--attributes_project', 'short': '-a'}}
  mosaicConfig   = mosaic_config.mosaicConfigFile(args.config)
  mosaicConfig   = mosaic_config.commandLineArguments(mosaicConfig, mosaicRequired)

  # Read the mosaicJson file to get information on how to process different annotations
  mosaicInfo = mosr.readMosaicJson(args.mosaic_json, args.reference)

  # Open an output tsv file to write annotations to
  outputFile = open(args.output_tsv, 'w')

  # Get the annotations from Mosaic in order to get the uids. OMIM annotations are private and so wont' be available
  # in the resource file. There are so many unique values that being a public annotation is a performance problem
  annotations       = api_va.getAnnotationDictNameIdUid(mosaicConfig, args.project_id)
  mosaicAnnotations = {}
  for annotationName in annotations: mosaicAnnotations[annotationName] = annotations[annotationName]['uid']

  # Loop over the annotations that are to be uploaded to Mosaic for this resource and get the annotation name, uid and
  # type
  annotations = {}
  uids        = []
  for annotation in mosaicInfo['resources']['OMIM']['annotations']:
    uid = mosaicAnnotations[annotation] if annotation in mosaicAnnotations else fail('Annotation ' + str(annotation) + ' is not present in the Mosaic project')
    tagType = mosaicInfo['resources']['OMIM']['annotations'][annotation]['type']
    annotations[annotation] = {'uid': uid, 'type': tagType}

  # Write the header line to the tsv file
  print('CHROM\tSTART\tEND\tREF\tALT\t', '\t'.join(str(x['uid']) for x in annotations.values()), sep = '', file = outputFile)

  # Loop over all records in the vcf file
  for record in os.popen(bcftools.query(bcftoolsExe, args.input_vcf, annotations.keys())).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')
  
    # If all values are '.', this line can be ignored
    uniqueValues = set(fields[5:])
    if len(uniqueValues) == 1 and list(uniqueValues)[0] == '.': continue

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

  # Close the output tsv file
  outputFile.close()

# Input options
def parseCommandLine():
  global version
  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  parser.add_argument('--input_vcf', '-i', required = True, metavar = 'string', help = 'The input vcf file to annotate')
  parser.add_argument('--output_tsv', '-o', required = True, metavar = 'string', help = 'The output tsv file')
  parser.add_argument('--reference', '-r', required = True, metavar = 'string', help = 'The genome reference file used')
  parser.add_argument('--tools_directory', '-s', required = True, metavar = 'string', help = 'The path to the directory where the tools live')
  parser.add_argument('--utils_directory', '-l', required = True, metavar = 'string', help = 'The path to the public-utils directory')
  parser.add_argument('--mosaic_json', '-m', required = True, metavar = 'string', help = 'The json file describing the Mosaic parameters')
  parser.add_argument('--project_id', '-p', required = True, metavar = 'string', help = 'The project id that variants will be uploaded to')

  # Mosaic arguments
  parser.add_argument('--config', '-c', required = True, metavar = "string", help = "The config file for Mosaic")
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
  print(message, sep = "")
  exit(1)

# Initialise global variables

# Pipeline version
version = "0.0.1"

# Mosaic aip information
mosaicConfig = {}

# Define the bcftools executable
bcftoolsExe = False

# Define the allowed annotation classes
allowedClasses = ['A', 'B', 'C', 'clinvar', 'compound', 'OMIM']

if __name__ == "__main__":
  main()
