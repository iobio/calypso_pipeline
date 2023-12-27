from os.path import exists

import argparse
import os
import json
import subprocess
import sys

# Add the path of the common functions and import them
from sys import path
path.append(os.path.dirname(__file__) + '/../components')
import calypso_path as cpath

def main():
  global annIds
  global mosaicConfig

  # Parse the command line
  args = parseCommandLine()

  # Build the path to include all components
  sys.path = cpath.buildPath(sys.path, args.utils_directory)
  import mosaic_config
  import api_variant_annotations as api_va
  import api_variant_filters as api_vf

  # Parse the Mosaic config file to get the token and url for the api calls
  mosaicRequired = {'MOSAIC_TOKEN': {'value': args.token, 'desc': 'An access token', 'long': '--token', 'short': '-t'},
                    'MOSAIC_URL': {'value': args.url, 'desc': 'The api url', 'long': '--url', 'short': '-u'}}
  mosaicConfig   = mosaic_config.mosaicConfigFile(args.config)
  mosaicConfig   = mosaic_config.commandLineArguments(mosaicConfig, mosaicRequired)

  # Variant filters will be added to the project and we need to define the annotation ids for the displayed columns. This will be
  # the Exomiser annotations just created as well as ClinVar and the gene symbol, the consequence and the genotypes. We need to
  # get the ids for some of these columns
  annotationIds = api_va.getAnnotationDictNameIdUid(mosaicConfig, args.project_id)

  # The following private variant annotations need to be created in Mosaic. These need to be private as they
  # are ranks or scores that are specific to the project they are run in - the same variant could have a 
  # different rank in a different project, and so these cannot be public, as only one value would be allowed
  # across projects:
  #
  # Exomiser rank
  # Exomiser p-value
  # Exomiser gene combined score
  # Exomiser gene phenotype score
  # Exomiser gene variant score
  # Exomiser variant score
  # Exomiser mode of inheritance
  createAnnotation(api_va, args.project_id, annotationIds, 'Exomiser Rank', 'rank', 'float')
  createAnnotation(api_va, args.project_id, annotationIds, 'Exomiser P-Value', 'pvalue', 'float')
  createAnnotation(api_va, args.project_id, annotationIds, 'Exomiser Gene Combined Score', 'comb', 'float')
  createAnnotation(api_va, args.project_id, annotationIds, 'Exomiser Gene Phenotype Score', 'pheno', 'float')
  createAnnotation(api_va, args.project_id, annotationIds, 'Exomiser Gene Variant Score', 'geneVar', 'float')
  createAnnotation(api_va, args.project_id, annotationIds, 'Exomiser Variant Score', 'variant', 'float')
  createAnnotation(api_va, args.project_id, annotationIds, 'Exomiser MOI', 'moi', 'string')

###############
###############
############### THIS WILL NEED UPDATING TO HANDLE DIFFERENT REFERENCES AND WHEN WE MOVE TO A NEW GENE SYMBOL
###############
###############
  displayColumnUids = ['"' + annotationIds['Gene Symbol GRCh38']['uid'] + '"', '"' + annotationIds['Gene Consequence GRCh38']['uid'] + '"', '"' + annotationIds['Genotype']['uid'] + '"']
  displayColumnUids.append('"' + annIds['rank']['uid'] + '"')
  displayColumnUids.append('"' + annIds['pvalue']['uid'] + '"')
  displayColumnUids.append('"' + annIds['comb']['uid'] + '"')
  displayColumnUids.append('"' + annIds['pheno']['uid'] + '"')
  displayColumnUids.append('"' + annIds['geneVar']['uid'] + '"')
  displayColumnUids.append('"' + annIds['variant']['uid'] + '"')
  displayColumnUids.append('"' + annIds['moi']['uid'] + '"')

  # Open the output file
  outputFile = open(args.output, 'w')
  print('CHROM', 'START', 'END', 'REF', 'ALT', annIds['rank']['uid'], annIds['pvalue']['uid'], annIds['comb']['uid'], annIds['pheno']['uid'], annIds['geneVar']['uid'], annIds['variant']['uid'], annIds['moi']['uid'], sep = '\t', file = outputFile)

  # Store the modes of inheritance that have been seen. These will be used to defined variant filters
  categories = {}
  categories['Top Candidates'] = []
  categories['All Candidates']   = []

  # Store the variants information
  variants    = {}
  variantInfo = {}

  # Check that the exomiser output file exists and open it for parsing
  if not exists(args.input): fail('The file ' + str(args.input) + ' was not found')
  variantsFile = open(args.input, 'r')
  for i, variant in enumerate(variantsFile.readlines()):

    # Skip the header line
    if variant.startswith('#'): continue

    # Split the line on tab and extract the values of interest
    fields = variant.rstrip().split('\t')
    chrom  = fields[14]
    start  = int(fields[15])
    end    = int(fields[16]) + 1
    rank   = fields[0]
    variantInfo[i] = {'chrom': chrom, 
                   'start': start,
                   'end': end,
                   'ref': fields[17],
                   'alt': fields[18],
                   'rank': rank,
                   'gene': fields[2],
                   'pvalue': fields[5],
                   'comb': fields[6],
                   'pheno': fields[7],
                   'geneVar': fields[8],
                   'variant': fields[9],
                   'moi': fields[4]}

    # If the ref or alt allele are 'N', change them to '*'
    if ref == 'N': ref = '*'
    if alt == 'N': alt = '*'

    # Multiple variants can have the same rank, and the same variant can have multiple ranks (e.g. based on different modes
    # of inheritance). Store the variant information so that these can be consolidated before being exported
    if chrom not in variants: variants[chrom] = {}
    if start not in variants[chrom]: variants[chrom][start] = {}
    if end not in variants[chrom][start]:
      variants[chrom][start][end] = {'ref': variantInfo[i]['ref'], 'alt': variantInfo[i]['alt'], 'rank': rank, 'pvalue': variantInfo[i]['pvalue'], 'comb': variantInfo[i]['comb'], 'pheno': variantInfo[i]['pheno'], 'geneVar': variantInfo[i]['geneVar'], 'variant': variantInfo[i]['variant'], 'moi': variantInfo[i]['moi']}
    else:
      variants[chrom][start][end]['rank']    = variants[chrom][start][end]['rank'] + ',' + rank
      variants[chrom][start][end]['pvalue']  = variants[chrom][start][end]['pvalue'] + ',' + variantInfo[i]['pvalue']
      variants[chrom][start][end]['comb']    = variants[chrom][start][end]['comb'] + ',' + variantInfo[i]['comb']
      variants[chrom][start][end]['pheno']   = variants[chrom][start][end]['pheno'] + ',' + variantInfo[i]['pheno']
      variants[chrom][start][end]['geneVar'] = variants[chrom][start][end]['geneVar'] + ',' + variantInfo[i]['geneVar']
      variants[chrom][start][end]['variant'] = variants[chrom][start][end]['variant'] + ',' + variantInfo[i]['variant']
      variants[chrom][start][end]['moi']     = variants[chrom][start][end]['moi'] + ',' + variantInfo[i]['moi']

    # If this is a new mode of inheritance, create a new category and store the variant id with it, unless
    # the p-value is less than the cut-off, in which case, the variant should go in the top candidates
    # category
    if variantInfo[i]['moi'] not in categories: categories[variantInfo[i]['moi']] = []
    if float(variantInfo[i]['pvalue']) < float(args.pvalue): categories['Top Candidates'].append(i)
    else:
      categories['All Candidates'].append(i)
      categories[variantInfo[i]['moi']].append(i)

  # Loop over all variants and write them to file
  for chrom in variants:
    for start in variants[chrom]:
      for end in variants[chrom][start]:
         print(chrom, start, end, variants[chrom][start][end]['ref'], variants[chrom][start][end]['alt'], variants[chrom][start][end]['rank'], variants[chrom][start][end]['pvalue'], variants[chrom][start][end]['comb'], variants[chrom][start][end]['pheno'], variants[chrom][start][end]['geneVar'], variants[chrom][start][end]['variant'], variants[chrom][start][end]['moi'], sep = '\t', file = outputFile)

  # Create an Exomiser category of variant filters and create filters for the top candidates and the different
  # modes of inheritance. First, get all the existing filters and store those in the Exomiser category
  existingFilters = {}
  data            = api_vf.getVariantFilters(mosaicConfig, args.project_id)
  for variantFilter in data:
    if variantFilter['category'] == 'Exomiser': existingFilters[variantFilter['name']] = variantFilter['id']

  # Loop over the variant categories (top candidates and modes of inheritance) and create the filters in the 'Exomiser' category if
  # they don't already exist, otherwise update them
  for name in categories:
    annotationFilters = defineAnnotationFilters(name, annIds, args.pvalue, name)
    if name in existingFilters: api_vf.updateVariantFilterColSort(mosaicConfig, args.project_id, name, existingFilters[name], displayColumnUids, annIds['rank']['uid'], 'ascending', annotationFilters)
    else: filterId = api_vf.createVariantFilterWithDisplay(mosaicConfig, args.project_id, name, 'Exomiser', displayColumnUids, annIds['rank']['uid'], 'ascending', annotationFilters)

  # Close the input and output files
  variantsFile.close()
  outputFile.close()

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
  parser.add_argument('--output', '-o', required = True, metavar = 'string', help = 'The output tsv file to upload to Mosaic')
  parser.add_argument('--pvalue', '-e', required = True, metavar = 'float', help = 'The p-value cutoff for variants to appear in the top candidates list')
  parser.add_argument('--project_id', '-p', required = True, metavar = 'float', help = 'The Mosaic project id that these exomiser results will be uploaded to')

  # The public-utils directory
  parser.add_argument('--utils_directory', '-l', required = True, metavar = 'string', help = 'The path to the public-utils directory')

  # Version
  parser.add_argument('--version', '-v', action='version', version='Exomiser variant output parser version: ' + str(version))

  return parser.parse_args()

# Check if the exomiser annotations exist and if not, create them
def createAnnotation(api_va, projectId, annotationIds, name, tag, annType):
  global annIds
  global mosaicConfig

  # If the annotation already exists, 
  if name in annotationIds:
    annId  = annotationIds[name]['id']
    annUid = annotationIds[name]['uid']

  # Otherwise, create the annotation
  else:
    data   = api_va.createPrivateAnnotationCategoryIdUid(mosaicConfig, name, annType, projectId, 'Exomiser')
    annId  = data['id']
    annUid = data['uid']

  # Update the annIds array with the annotation id and uid
  annIds[tag] = {'id': annId, 'uid': annUid}

# Define the annotation filters
def defineAnnotationFilters(filterName, annotationIds, pValue, name):

  # Handle the Top Candidates filter
  if str(filterName) == 'Top Candidates':
    jsonFilters = {
      "annotation_filters": [
        {
          "uid": annotationIds['pvalue']['uid'],
          "max": pValue,
          "include_nulls": False
        }
      ]
    }

  # The remaining filters are all based on the mode of inheritance
  elif str(filterName) == 'All Candidates':
    jsonFilters = {
      "annotation_filters": [
        {
          "uid": annotationIds['rank']['uid'],
          "min": "1",
          "include_nulls": False
        }
      ]
    }

  # The remaining filters are all based on the mode of inheritance
  else:
    jsonFilters = {
      "annotation_filters": [
        {        
          "uid": annotationIds['rank']['uid'],
          "min": "1",
          "include_nulls": False
        },
        {        
          "uid": annotationIds['moi']['uid'],
          "values": [name],
          "include_nulls": False
        }
      ]
    }

  # Return the json object
  return jsonFilters

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
