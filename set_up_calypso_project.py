from __future__ import print_function

import os
import argparse
import json
import math
import glob
import importlib
import sys

from datetime import date
from os.path import exists
from sys import path

# Add the path of the common functions and import them
path.append(os.path.dirname(__file__) + '/components')
import calypso_annotations as anno
import calypso_bash_script as bashScript
import clinvar_compare as cv
import calypso_mosaic as mos
import calypso_mosaic_resources as mosr
import calypso_path as cpath
import calypso_reanalysis as rean
import calypso_resources as res
import calypso_samples as sam
import calypso_slivar_trio as slivarTrio
import calypso_toml as toml
import calypso_upload as upload
import calypso_vcf_files as vcfs
import tools_bcftools as bcftools

def main():
  global mosaicConfig
  global workingDir
  global version
  global allowedReferences

  print()
  print('Starting the Calypso pipeline')

  # Parse the command line
  args = parseCommandLine()

  # Check the supplied parameters are as expected, then expand the path to enable scripts and api commands
  # to be accessed for Calypso
  checkArguments(args)
  sys.path = cpath.buildPath(sys.path, args.utils_directory)
  rootPath = os.path.dirname(__file__)

  # Import additional tools
  import mosaic_config
  import api_pedigree as api_ped
  import api_project_attributes as api_pa
  import api_project_settings as api_ps
  import api_samples as api_s
  import api_sample_attributes as api_sa
  import api_sample_files as api_sf
  import api_sample_hpo_terms as api_shpo
  import api_variants as api_v
  import api_variant_annotations as api_va
  import api_variant_filters as api_vf
  import add_variant_filters as pu_avf
  import variant_filters as vFilters

  # Parse the Mosaic config file to get the token and url for the api calls
  mosaicRequired = {'MOSAIC_TOKEN': {'value': args.token, 'desc': 'An access token', 'long': '--token', 'short': '-t'},
                    'MOSAIC_URL': {'value': args.url, 'desc': 'The api url', 'long': '--url', 'short': '-u'},
                    'MOSAIC_ATTRIBUTES_PROJECT_ID': {'value': args.attributes_project, 'desc': 'The public attribtes project id', 'long': '--attributes_project', 'short': '-a'}}
  mosaicConfig   = mosaic_config.mosaicConfigFile(args.config)
  mosaicConfig   = mosaic_config.commandLineArguments(mosaicConfig, mosaicRequired)

  # Get the project reference
  reference = api_ps.getProjectReference(mosaicConfig, args.project_id) if not args.reference else args.reference
  if reference not in allowedReferences: fail('The specified reference (' + str(reference) + ') is not recognised. Allowed values are: ' + str(', '.join(allowedReferences)))
  print('  Using the reference:                  ', reference, sep = '')

  # Read the resources json file to identify all the resources that will be used in the annotation pipeline
  #if args.resource_json: resourceInfo = res.checkResources(reference, args.data_directory, args.tools_directory, args.resource_json)
  #else: resourceInfo = res.checkResources(reference, args.data_directory, args.tools_directory, args.data_directory + 'resources_' + str(reference) + '.json')
  #resourceInfo = res.readResources(reference, rootPath, resourceInfo, args.no_vep)

  # Define the tools to be used by Calypso
  #resourceInfo = res.calypsoTools(resourceInfo)

  # Define paths to be used by Calypso
  #setWorkingDir(resourceInfo['version'])

  # Get all the project attributes in the target Mosaic project
  #projectAttributes = mos.getProjectAttributes(mosaicConfig, args.project_id, api_pa)
  #publicAttributes  = mos.getPublicProjectAttributes(mosaicConfig, api_pa)

  # Read the Mosaic json and validate its contents
  if not args.mosaic_json: args.mosaic_json = args.data_directory + 'resources_mosaic_' + reference + '.json'
  mosaicInfo         = mosr.readMosaicJson(args.mosaic_json, reference)
  projectAnnotations = mos.getProjectAnnotations(mosaicConfig, args.project_id, api_va)
  publicAnnotations  = mos.getPublicAnnotations(mosaicConfig, args.project_id, api_va)
  privateAnnotations, projectAnnotations = mos.createPrivateAnnotations(mosaicConfig, mosaicInfo['resources'], projectAnnotations, samples, args.project_id, api_va)
  mosr.checkPublicAnnotations(mosaicConfig, mosaicInfo['resources'], publicAnnotations, args.project_id, api_va)

  # Prepare the target project. This includes:
  # 1. Remove any unnecessary annotations from the project
  # 2. Import public annotations into the project
  # 3. Set the default annotations for the project
  # 4. Set up required variant filters
  # 5. Add and set project attributes
  if 'remove' in mosaicInfo: mos.removeAnnotations(mosaicConfig, mosaicInfo['remove'], projectAnnotations, args.project_id, api_va)
  projectAnnotations = mos.importAnnotations(mosaicConfig, mosaicInfo['resources'], projectAnnotations, publicAnnotations, args.project_id, api_va)
  mos.defaultAnnotations(mosaicConfig, mosaicInfo['defaultAnnotations'], publicAnnotations, privateAnnotations, args.project_id, api_ps)

  # Determine all of the variant filters (from the calypso_mosaic_filters.json) that are to be added; remove any filters that already
  # exist with the same name; fill out variant filter details not in the json (e.g. the uids of private annotations created by
  # Calypso); create the filters; and finally update the project settings to put the filters in the correct category and sort order.
  # Note that the filters to be applied depend on the family structure. E.g. de novo filters won't be added to projects without parents
  vFilters.setVariantFilters(mosaicConfig, api_ps, api_va, api_vf, args.project_id, args.variant_filters, mosaicSamples)

# Input options
def parseCommandLine():
  global version
  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  parser.add_argument('--data_directory', '-d', required = True, metavar = 'string', help = 'The path to the directory where the resources live')
  parser.add_argument('--tools_directory', '-s', required = False, metavar = 'string', help = 'The path to the directory where the tools to use live')
  parser.add_argument('--utils_directory', '-l', required = True, metavar = 'string', help = 'The path to the public-utils directory')
  parser.add_argument('--input_vcf', '-i', required = False, metavar = 'string', help = 'The input vcf file to annotate')
  parser.add_argument('--ped', '-e', required = False, metavar = 'string', help = 'The pedigree file for the family. Not required for singletons')
  parser.add_argument('--project_id', '-p', required = True, metavar = 'string', help = 'The project id that variants will be uploaded to')
  parser.add_argument('--config', '-c', required = True, metavar = 'string', help = 'The config file for Mosaic')
  parser.add_argument('--variant_filters', '-f', required = True, metavar = 'string', help = 'The json file describing the variant filters to apply to each project')

  # Optional pipeline arguments
  parser.add_argument('--reference', '-r', required = False, metavar = 'string', help = 'The reference genome to use. Allowed values: ' + ', '.join(allowedReferences))
  parser.add_argument('--resource_json', '-j', required = False, metavar = 'string', help = 'The json file describing the annotation resources')
  parser.add_argument('--no_vep', '-n', required = False, action = "store_false", help = "If set, do not run VEP annotation")
#  parser.add_argument('--no_comp_het', '-n', required = False, action = "store_true", help = "If set, comp het determination will be skipped")

  # Optional argument to handle HPO terms
  parser.add_argument('--hpo', '-o', required = False, metavar = "string", help = "A comma separate list of hpo ids for the proband")


  # Optional mosaic arguments
  parser.add_argument('--token', '-t', required = False, metavar = "string", help = "The Mosaic authorization token")
  parser.add_argument('--url', '-u', required = False, metavar = "string", help = "The base url for Mosaic")
  parser.add_argument('--attributes_project', '-a', required = False, metavar = "integer", help = "The Mosaic project id that contains public attributes")
  parser.add_argument('--mosaic_json', '-m', required = False, metavar = 'string', help = 'The json file describing the Mosaic parameters')

  # Version
  parser.add_argument('--version', '-v', action="version", version='Calypso annotation pipeline version: ' + str(version))

  return parser.parse_args()

# Check the supplied arguments
def checkArguments(args):

  # Ensure the data path ends with a "/", then add the reference directory
  if args.data_directory[-1] != "/": args.data_directory += "/"

  # Check that the public-utils directory exists and add to the path so scripts from here can be used
  if args.utils_directory[-1] != "/": args.utils_directory += "/"
  if not os.path.exists(args.utils_directory): fail('public-utils directory could not be found')

# Create a directory where all Calypso associated files will be stored
def setWorkingDir(version):
  global workingDir

  workingDir += version + "/"

  # If the directory doesn't exist, create it
  if not os.path.exists(workingDir): os.makedirs(workingDir)

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)

# Initialise global variables

# Pipeline version
version = "1.1.5"
date    = str(date.today())

# The working directory where all created files are kept
workingDir = os.getcwd() + "/calypso_v" + version + "r"

# Store information related to Mosaic
mosaicConfig = {}

# Store the allowed references that can be specified on the command line
allowedReferences = ['GRCh37', 'GRCh38']

if __name__ == "__main__":
  main()