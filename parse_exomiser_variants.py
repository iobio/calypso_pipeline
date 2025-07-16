from os.path import exists
from pprint import pprint
from datetime import datetime

import argparse
import os
import json
import subprocess
import sys

from sys import path

def main():

  # Parse the command line
  args = parseCommandLine()

  # Import the api client
  path.append(args.api_client)
  from mosaic import Mosaic, Project, Store
  api_store  = Store(config_file = args.client_config)
  api_mosaic = Mosaic(config_file = args.client_config)

  # Open an api client project object for the defined project
  project = api_mosaic.get_project(args.project_id)

  # Get the date
  date = datetime.today().strftime('%Y-%m-%d')

  # Variant filters will be added to the project and we need to define the annotation ids for the displayed columns. This will be
  # the Exomiser annotations just created as well as ClinVar and the gene symbol, the consequence and the genotypes. We need to
  # get the ids for some of these columns
  annotations = {}
  for annotation in project.get_variant_annotations():
    annotations[annotation['name']] = {'id': annotation['id'], 'uid': annotation['uid'], 'version_ids': {}}
    for version in annotation['annotation_versions']:
      annotations[annotation['name']]['version_ids'][version['version']] = version['id']

  # The following private variant annotations need to be created in Mosaic. These need to be private as they
  # are ranks or scores that are specific to the project they are run in - the same variant could have a 
  # different rank in a different project, and so these cannot be public, as only one value would be allowed
  # across projects:
  #
  # Exomiser rank
  # Exomiser contributing variant
  # Exomiser p-value
  # Exomiser gene combined score
  # Exomiser gene phenotype score
  # Exomiser gene variant score
  # Exomiser variant score
  # Exomiser mode of inheritance
  version_name = args.version_name if args.version_name else False
  annotations = create_annotation(project, args.project_id, date, annotations, 'Exomiser Rank', 'float', version_name, args.overwrite_versions)
  annotations = create_annotation(project, args.project_id, date, annotations, 'Exomiser Contributing Variant', 'float', version_name, args.overwrite_versions)
  annotations = create_annotation(project, args.project_id, date, annotations, 'Exomiser P-Value', 'float', version_name, args.overwrite_versions)
  annotations = create_annotation(project, args.project_id, date, annotations, 'Exomiser Gene Combined Score', 'float', version_name, args.overwrite_versions)
  annotations = create_annotation(project, args.project_id, date, annotations, 'Exomiser Gene Phenotype Score', 'float', version_name, args.overwrite_versions)
  annotations = create_annotation(project, args.project_id, date, annotations, 'Exomiser Gene Variant Score', 'float', version_name, args.overwrite_versions)
  annotations = create_annotation(project, args.project_id, date, annotations, 'Exomiser Variant Score', 'float', version_name, args.overwrite_versions)
  annotations = create_annotation(project, args.project_id, date, annotations, 'Exomiser MOI', 'string', version_name, args.overwrite_versions)

  # Open the output file
  output_file = open(args.output, 'w')
  print('CHROM', 'START', 'END', 'REF', 'ALT', \
        str(annotations['Exomiser Rank']['uid']) + '@' + str(annotations['Exomiser Rank']['version']), \
        str(annotations['Exomiser P-Value']['uid']) + '@' + str(annotations['Exomiser P-Value']['version']), \
        str(annotations['Exomiser Gene Combined Score']['uid']) + '@' + str(annotations['Exomiser Gene Combined Score']['version']), \
        str(annotations['Exomiser Gene Phenotype Score']['uid']) + '@' + str(annotations['Exomiser Gene Phenotype Score']['version']), \
        str(annotations['Exomiser Gene Variant Score']['uid']) + '@' + str(annotations['Exomiser Gene Variant Score']['version']), \
        str(annotations['Exomiser Variant Score']['uid']) + '@' + str(annotations['Exomiser Variant Score']['version']), \
        str(annotations['Exomiser MOI']['uid']) + '@' + str(annotations['Exomiser MOI']['version']), \
        str(annotations['Exomiser Contributing Variant']['uid']) + '@' + str(annotations['Exomiser Contributing Variant']['version']), \
        sep = '\t', file = output_file)

  # Store the variants information
  variants = {}
  variant_info = {}

  # Check that the exomiser output file exists and open it for parsing
  if not exists(args.input):
    fail('The file ' + str(args.input) + ' was not found')
  variants_file = open(args.input, 'r')
  for i, variant in enumerate(variants_file.readlines()):

    # Skip the header line
    if variant.startswith('#'):
      continue

    # Split the line on tab and extract the values of interest
    fields = variant.rstrip().split('\t')
    chrom  = fields[14]
    start  = int(fields[15])
    end    = int(fields[16]) + 1
    rank   = fields[0]
    variant_info[i] = {'chrom': chrom, 
                   'start': start,
                   'end': end,
                   'ref': fields[17],
                   'alt': fields[18],
                   'Exomiser Rank': rank,
                   'gene': fields[2],
                   'Exomiser P-Value': fields[5],
                   'Exomiser Gene Combined Score': fields[6],
                   'Exomiser Gene Phenotype Score': fields[7],
                   'Exomiser Gene Variant Score': fields[8],
                   'Exomiser Variant Score': fields[9],
                   'Exomiser MOI': fields[4],
                   'Exomiser Contributing Variant': fields[10]}

    # If the ref or alt allele are 'N', change them to '*'
    if variant_info[i]['ref'] == 'N':
      variant_info[i]['ref'] = '*'
    if variant_info[i]['alt'] == 'N':
      variant_info[i]['alt'] = '*'

    # Multiple variants can have the same rank, and the same variant can have multiple ranks (e.g. based on different modes
    # of inheritance). Store the variant information so that these can be consolidated before being exported
    if chrom not in variants:
      variants[chrom] = {}
    if start not in variants[chrom]:
      variants[chrom][start] = {}
    if end not in variants[chrom][start]:
      variants[chrom][start][end] = {'ref': variant_info[i]['ref'], \
                                     'alt': variant_info[i]['alt'], \
                                     'Exomiser Rank': rank,
                                     'Exomiser P-Value': variant_info[i]['Exomiser P-Value'], \
                                     'Exomiser Gene Combined Score': variant_info[i]['Exomiser Gene Combined Score'], \
                                     'Exomiser Gene Phenotype Score': variant_info[i]['Exomiser Gene Phenotype Score'], \
                                     'Exomiser Gene Variant Score': variant_info[i]['Exomiser Gene Variant Score'], \
                                     'Exomiser Variant Score': variant_info[i]['Exomiser Variant Score'], \
                                     'Exomiser MOI': variant_info[i]['Exomiser MOI'], \
                                     'Exomiser Contributing Variant': variant_info[i]['Exomiser Contributing Variant']}
    else:
      variants[chrom][start][end]['Exomiser Rank'] = variants[chrom][start][end]['Exomiser Rank'] + ',' + variant_info[i]['Exomiser Rank']
      variants[chrom][start][end]['Exomiser P-Value'] = variants[chrom][start][end]['Exomiser P-Value'] + ',' + variant_info[i]['Exomiser P-Value']
      variants[chrom][start][end]['Exomiser Gene Combined Score'] = variants[chrom][start][end]['Exomiser Gene Combined Score'] + ',' + variant_info[i]['Exomiser Gene Combined Score']
      variants[chrom][start][end]['Exomiser Gene Phenotype Score'] = variants[chrom][start][end]['Exomiser Gene Phenotype Score'] + ',' + variant_info[i]['Exomiser Gene Phenotype Score']
      variants[chrom][start][end]['Exomiser Gene Variant Score'] = variants[chrom][start][end]['Exomiser Gene Variant Score'] + ',' + variant_info[i]['Exomiser Gene Variant Score']
      variants[chrom][start][end]['Exomiser Variant Score'] = variants[chrom][start][end]['Exomiser Variant Score'] + ',' + variant_info[i]['Exomiser Variant Score']
      variants[chrom][start][end]['Exomiser MOI'] = variants[chrom][start][end]['Exomiser MOI'] + ',' + variant_info[i]['Exomiser MOI']
      variants[chrom][start][end]['Exomiser Contributing Variant'] = variants[chrom][start][end]['Exomiser Contributing Variant'] + ',' + variant_info[i]['Exomiser Contributing Variant']

  # Loop over all variants and write them to file
  for chrom in variants:
    for start in variants[chrom]:
      for end in variants[chrom][start]:
        print(chrom, start, end, 
              variants[chrom][start][end]['ref'], \
              variants[chrom][start][end]['alt'], \
              variants[chrom][start][end]['Exomiser Rank'], \
              variants[chrom][start][end]['Exomiser P-Value'], \
              variants[chrom][start][end]['Exomiser Gene Combined Score'], \
              variants[chrom][start][end]['Exomiser Gene Phenotype Score'], \
              variants[chrom][start][end]['Exomiser Gene Variant Score'], \
              variants[chrom][start][end]['Exomiser Variant Score'], \
              variants[chrom][start][end]['Exomiser MOI'], \
              variants[chrom][start][end]['Exomiser Contributing Variant'], sep = '\t', file = output_file)

  # Close the input and output files
  variants_file.close()
  output_file.close()

# Input options
def parseCommandLine():
  global version

  parser = argparse.ArgumentParser(description='Process the command line')

  # Define the location of the api_client and the ini config file
  parser.add_argument('--api_client', '-a', required = True, metavar = 'string', help = 'The api_client directory')
  parser.add_argument('--client_config', '-c', required = True, metavar = 'string', help = 'The ini config file for Mosaic')

  # Required arguments
  parser.add_argument('--input', '-i', required = True, metavar = 'string', help = 'The exomiser variants.tsv file')
  parser.add_argument('--output', '-o', required = True, metavar = 'string', help = 'The output tsv file to upload to Mosaic')
  parser.add_argument('--project_id', '-p', required = True, metavar = 'float', help = 'The Mosaic project id that these exomiser results will be uploaded to')

  # Optional arguments
  parser.add_argument('--version_name', '-v', required = False, metavar = 'string', help = 'Optionally supply the name of the annotation version')
  parser.add_argument('--overwrite_versions', '-ov', required = False, action = 'store_true', help = 'If set, existing annotation versions with the same name will be overwritten')

  return parser.parse_args()

# Check if the exomiser annotations exist and if not, create them
def create_annotation(project, project_id, date, annotations, name, annotation_type, version_name, overwrite):
  version_id = False

  # If the annotation doesn't exist, create it. If it does exist and the option to keep existing annotations
  # has been selected create a new annotation, that includes the date
  if not name in annotations:
    version_name = str(date) if not version_name else version_name
    data = project.post_variant_annotation(name = name, category = 'Exomiser', value_type = annotation_type, privacy_level = 'private', value_truncate_type = 'middle')
    annotation_id = data['id']
    annotation_uid = data['uid']
    version_id = project.post_create_annotation_version(annotation_id, version_name)['id']

    # Update the annotations array with the annotation id and uid
    annotations[name] = {'id': annotation_id, 'uid': annotation_uid, 'version_id': version_id, 'version': version_name}

  # If the annotation exists, create a new version
  else:
    version_name = str(date) if not version_name else version_name
    if version_name in annotations[name]['version_ids']:

      # If the version can be overwritten, delete the existing annotation version, so a new one with this
      # name can be created
      if overwrite:
        try:
          project.delete_variant_annotation_version(annotations[name]['id'], annotations[name]['version_ids'][version_name])
        except Exception as e:
          fail('Couldn\'t delete annotation version in. Error is: ' + str(e))

      # Otherwise fail, as the version already exists
      else:
        fail('A version with today\'s date already exists for annotation with the name ' + str(name))
    version_id = project.post_create_annotation_version(annotations[name]['id'], version_name)['id']

    # Update the annotations array with the annotation id and uid
    annotations[name]['version_ids']['version_name'] = version_id
    annotations[name]['version'] = version_name

    # Set the latest version of the annotation to the version just created
    project.put_variant_annotation(annotations[name]['id'], name = name, latest_version_id = version_id)
 
  # Return the updated annotations
  return annotations

# If the script fails, provide an error message and exit
def fail(message):
  print('ERROR: ', message, sep = '')
  exit(1)

# Initialise global variables

if __name__ == "__main__":
  main()
