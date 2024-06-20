from os.path import exists

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

  # Variant filters will be added to the project and we need to define the annotation ids for the displayed columns. This will be
  # the Exomiser annotations just created as well as ClinVar and the gene symbol, the consequence and the genotypes. We need to
  # get the ids for some of these columns
  #
  # For now, just use the default annotation version
  annotations = {}
  for annotation in project.get_variant_annotations():
    for version in annotation['annotation_versions']:
      if version['version'] == 'default':
        version_id = version['id']
        break
    annotations[annotation['name']] = {'id': annotation['id'], 'uid': annotation['uid'], 'version_id': version_id}

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
  annotations = create_annotation(project, args.project_id, annotations, 'Exomiser Rank', 'float')
  annotations = create_annotation(project, args.project_id, annotations, 'Exomiser Contributing Variant', 'float')
  annotations = create_annotation(project, args.project_id, annotations, 'Exomiser P-Value', 'float')
  annotations = create_annotation(project, args.project_id, annotations, 'Exomiser Gene Combined Score', 'float')
  annotations = create_annotation(project, args.project_id, annotations, 'Exomiser Gene Phenotype Score', 'float')
  annotations = create_annotation(project, args.project_id, annotations, 'Exomiser Gene Variant Score', 'float')
  annotations = create_annotation(project, args.project_id, annotations, 'Exomiser Variant Score', 'float')
  annotations = create_annotation(project, args.project_id, annotations, 'Exomiser MOI', 'string')

  # Get the project reference
  reference = project.get_project_settings()['reference']

  if reference == 'GRCh37':
    display_column_ids = [str(annotations['Gene Name']['version_id']), str(annotations['Gene Consequence GRCh37']['version_id']), str(annotations['Genotype']['version_id'])]
  elif reference == 'GRCh38':
    display_column_ids = [str(annotations['Gene Name']['version_id']), str(annotations['Gene Consequence GRCh38']['version_id']), str(annotations['Genotype']['version_id'])]
  else:
    fail('Project has an unknown reference: ' + reference)

  display_column_ids.append(str(annotations['Exomiser Rank']['version_id']))
  display_column_ids.append(str(annotations['Exomiser P-Value']['version_id']))
  display_column_ids.append(str(annotations['Exomiser Contributing Variant']['version_id']))
  display_column_ids.append(str(annotations['Exomiser Gene Combined Score']['version_id']))
  display_column_ids.append(str(annotations['Exomiser Gene Phenotype Score']['version_id']))
  display_column_ids.append(str(annotations['Exomiser Gene Variant Score']['version_id']))
  display_column_ids.append(str(annotations['Exomiser Variant Score']['version_id']))
  display_column_ids.append(str(annotations['Exomiser MOI']['version_id']))

  # Open the output file
  output_file = open(args.output, 'w')
  print('CHROM', 'START', 'END', 'REF', 'ALT', \
        annotations['Exomiser Rank']['uid'], \
        annotations['Exomiser P-Value']['uid'], \
        annotations['Exomiser Gene Combined Score']['uid'], \
        annotations['Exomiser Gene Phenotype Score']['uid'], \
        annotations['Exomiser Gene Variant Score']['uid'], \
        annotations['Exomiser Variant Score']['uid'], \
        annotations['Exomiser MOI']['uid'], \
        annotations['Exomiser Contributing Variant']['uid'], sep = '\t', file = output_file)

  # Store the modes of inheritance that have been seen. These will be used to defined variant filters
  categories = {}
  categories['Top Candidates'] = []
  categories['All Candidates']   = []

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

    # If this is a new mode of inheritance, create a new category and store the variant id with it, unless
    # the p-value is less than the cut-off, in which case, the variant should go in the top candidates
    # category
    if float(variant_info[i]['Exomiser P-Value']) < float(args.pvalue):
      categories['Top Candidates'].append(i)
    else:
      categories['All Candidates'].append(i)

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

  # Create an Exomiser category of variant filters and create filters for the top candidates and the different
  # modes of inheritance. First, get all the existing filters and store those in the Exomiser category
  existing_filters = {}
  for variant_filter in project.get_variant_filters():
    if variant_filter['category'] == 'Exomiser':
      existing_filters[variant_filter['name']] = variant_filter['id']

  # Loop over the variant categories (top candidates and modes of inheritance) and create the filters in the 'Exomiser' category if
  # they don't already exist, otherwise update them
  for name in categories:
    annotation_filters = define_annotation_filters(annotations, name, args.pvalue)
    if name in existing_filters:
      filter_id = existing_filters[name]
      project.update_variant_filter(filter_id, annotation_filters)
    else:
      sort_column_id = annotations['Exomiser Rank']['version_id']
      project.post_variant_filter(name = name, category = 'Exomiser', column_ids = display_column_ids, sort_column_id = sort_column_id, sort_direction = 'ascending', filter_data = annotation_filters)

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
  parser.add_argument('--pvalue', '-e', required = True, metavar = 'float', help = 'The p-value cutoff for variants to appear in the top candidates list')
  parser.add_argument('--project_id', '-p', required = True, metavar = 'float', help = 'The Mosaic project id that these exomiser results will be uploaded to')

  return parser.parse_args()

# Check if the exomiser annotations exist and if not, create them
def create_annotation(project, project_id, annotations, name, annotation_type):
  if not name in annotations:
    data = project.post_variant_annotation(name = name, category = 'Exomiser', value_type = annotation_type, privacy_level = 'private', value_truncate_type = 'middle')
    annotation_id = data['id']
    annotation_uid = data['uid']

    # Get the annotation version if of the default annotation
    for version in data['annotation_versions']:
      if version['version'] == 'default':
        version_id = version['id']
        break

    # Update the annotations array with the annotation id and uid
    annotations[name] = {'id': annotation_id, 'uid': annotation_uid, 'version_id': version_id}
 
  # Return the updated annotations
  return annotations

# Define the annotation filters
def define_annotation_filters(annotations, name, pvalue):

  # Handle the Top Candidates filter
  if str(name) == 'Top Candidates':
    json_filters = {
      "annotation_filters": [
        {
          "uid": annotations['Exomiser P-Value']['uid'],
          "annotation_version_id": annotations['Exomiser P-Value']['version_id'],
          "max": pvalue,
          "include_nulls": False
        }
      ]
    }

  # The remaining filters are all based on the mode of inheritance
  elif str(name) == 'All Candidates':
    json_filters = {
      "annotation_filters": [
        {
          "uid": annotations['Exomiser Rank']['uid'],
          "annotation_version_id": annotations['Exomiser Rank']['version_id'],
          "min": "1",
          "include_nulls": False
        }
      ]
    }

  # Return the json object
  return json_filters

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = '')
  exit(1)

# Initialise global variables

if __name__ == "__main__":
  main()
