#!/usr/bin/python

from __future__ import print_function
from datetime import date
from os.path import exists
from sys import path

import os
import argparse
import json
import math
import glob
import importlib

# Add the path of the common functions and import them
path.append(os.path.dirname(__file__) + '/components')
import calypso_annotations as anno
import calypso_bash_script as bashScript
import calypso_mosaic as mos
import calypso_mosaic_resources as mosr
import calypso_resources as res
import calypso_samples as sam
import calypso_toml as toml
import calypso_upload as upload
import calypso_variant_filters as vfilt

def main():
  global mosaicConfig
  global workingDir
  global version

  # Parse the command line
  args = parseCommandLine()

  # Check the supplied parameters are as expected, then expand the path to enable scripts and api commands
  # to be accessed for Calypso
  resourceInfo = checkArguments(args)
  buildPath(args)

  # Read the resources json file to identify all the resources that will be used in the annotation pipeline
  resourceInfo = res.readResources(args, resourceInfo)
  rootPath = os.path.dirname( __file__)
  setWorkingDir(resourceInfo['version'])

  # Import additional tools
  import mosaic_config
  import api_project_attributes as api_pa
  import api_project_settings as api_ps
  import api_samples as api_s
  import api_variants as api_v
  import api_variant_annotations as api_va
  import api_variant_filters as api_vf

  # Parse the Mosaic config file to get the token and url for the api calls
  mosaicRequired = {"token": True, "url": True, "attributesProjectId": True}
  mosaicConfig   = mosaic_config.parseConfig(args, mosaicRequired)

  # Determine the id of the proband, and the order of the samples in the vcf header
  proband, samples = sam.getProband(mosaicConfig, args.ped, args.family_type, args.project_id, api_s)
  sam.getSampleOrder(samples, args.input_vcf)

  # Get all the project attributes in the target Mosaic project
  projectAttributes = mos.getProjectAttributes(mosaicConfig, args.project_id, api_pa)
  publicAttributes  = mos.getPublicProjectAttributes(mosaicConfig, api_pa)

  # Read the Mosaic json and validate its contents
  if not args.mosaic_json: args.mosaic_json = args.data_directory + 'resources_mosaic_GRCh' + args.reference + '.json'
  mosaicInfo         = mosr.readMosaicJson(args.mosaic_json, args.reference)
  projectAnnotations = mos.getProjectAnnotations(mosaicConfig, args.project_id, api_va)
  publicAnnotations  = mos.getPublicAnnotations(mosaicConfig, args.project_id, api_va)
  privateAnnotations = mos.createPrivateAnnotations(mosaicConfig, mosaicInfo['resources'], projectAnnotations, samples, args.project_id, api_va)
  mosr.checkPublicAnnotations(mosaicConfig, mosaicInfo['resources'], publicAnnotations, args.project_id, api_va)

  # Build the toml file for vcfanno. This defines all of the annotations that vcfanno needs to use
  # with the path to the required files
  tomlFilename = toml.buildToml(workingDir, resourceInfo)

  # Generate the bash script to run the annotation pipeline
  bashFilename, bashFile = bashScript.openBashScript(workingDir)
  filteredVcf = bashScript.bashResources(resourceInfo, workingDir, bashFile, args.input_vcf, args.ped, tomlFilename)
  bashScript.samplesFile(bashFile)
  bashScript.cleanedVcf(bashFile)
  bashScript.annotateVcf(bashFile, resourceInfo)
  bashScript.filterVariants(bashFile, proband, resourceInfo)
  bashScript.rareDiseaseVariants(bashFile)
  bashScript.deleteFiles(bashFile)

  # Process the filtered vcf file to extract the annotations to be uploaded to Mosaic
  print('# Generate tsv files to upload annotations to Mosaic', file = bashFile)
  uploadFilename, uploadFile = upload.openUploadAnnotationsFile(workingDir)
  for resource in mosaicInfo['resources']:
    tsv = anno.createAnnotationTsv(mosaicInfo, resource, os.path.dirname(__file__) + '/components', args.reference, args.config, args.mosaic_json, bashFile)
    upload.uploadAnnotations(args.utils_directory, tsv, mosaicInfo['resources'][resource]['project_id'], args.config, uploadFile)
  upload.closeUploadAnnotationsFile(uploadFilename, uploadFile)

  # Close the bash script
  bashScript.finishScript(bashFile, bashFilename)

  # Prepare the target project. This includes:
  # 1. Remove any unnecessary annotations from the project
  # 2. Import public annotations into the project
  # 3. Set the default annotations for the project
  # 4. Set up required variant filters
  # 5. Add and set project attributes
  if 'remove' in mosaicInfo: mos.removeAnnotations(mosaicConfig, mosaicInfo['remove'], projectAnnotations, args.project_id, api_va)
  mos.importAnnotations(mosaicConfig, mosaicInfo['resources'], projectAnnotations, publicAnnotations, args.project_id, api_va)
  mos.defaultAnnotations(mosaicConfig, mosaicInfo['defaultAnnotations'], publicAnnotations, privateAnnotations, args.project_id, api_ps)
  vfilt.variantFilters(mosaicConfig, rootPath, samples, privateAnnotations, args.project_id, api_vf)
  mos.updateCalypsoAttributes(mosaicConfig, projectAttributes, publicAttributes, version, args.project_id, api_pa)

  # Generate scripts to upload filtered variants and annotations to Mosaic
  upload.uploadVariants(mosaicInfo, mosaicConfig, workingDir, args.config, filteredVcf, args.project_id, api_v)
  #upload.uploadAnnotations(args)

#  # If a previous Calypso vcf file was set, get the version assigned to the vcf
#  if args.previous_vcf: getPreviousVcfInfo(args)

  # Output a summary file listing the actions undertaken by Calypso with all version histories
  res.calypsoSummary(workingDir, version, resourceInfo, args.reference)

# Input options
def parseCommandLine():
  global version
  global allowedReferences
  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  parser.add_argument('--reference', '-r', required = True, metavar = 'string', help = 'The reference genome to use. Allowed values: ' + ', '.join(allowedReferences))
  parser.add_argument('--family_type', '-f', required = True, metavar = 'string', help = 'The familty structure. Allowed values: ' + ', '. join(allowedFamily))
  parser.add_argument('--data_directory', '-d', required = True, metavar = 'string', help = 'The path to the directory where the resources live')
  parser.add_argument('--utils_directory', '-l', required = True, metavar = 'string', help = 'The path to the public-utils directory')
  parser.add_argument('--input_vcf', '-i', required = True, metavar = 'string', help = 'The input vcf file to annotate')
  parser.add_argument('--ped', '-e', required = True, metavar = 'string', help = 'The pedigree file for the family. Not required for singletons')
  parser.add_argument('--project_id', '-p', required = True, metavar = 'string', help = 'The project id that variants will be uploaded to')
  parser.add_argument('--config', '-c', required = True, metavar = "string", help = "The config file for Mosaic")

  # Optional reanalysis arguments
#  parser.add_argument('--previous_vcf', '-s', required = False, metavar = "string", help = "A previously annotated Calypso vcf file")

  # Optional pipeline arguments
  parser.add_argument('--resource_json', '-j', required = False, metavar = 'string', help = 'The json file describing the annotation resources')
#  parser.add_argument('--no_comp_het', '-n', required = False, action = "store_true", help = "If set, comp het determination will be skipped")
#
  # Optional argument to handle HPO terms
#  parser.add_argument('--hpo', '-o', required = False, metavar = "string", help = "A comma separate list of hpo ids for the proband")

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
  global allowedReferences
  global allowedFamily
  resourceInfo = {}

  # Check the reference is allowed
  if args.reference not in allowedReferences:
    message = 'The allowed references (--reference, -r) are:'
    for ref in allowedReferences: message += '\n  ' + str(ref)
    fail(message)

  # Ensure the data path ends with a "/", then add the reference directory
  if args.data_directory[-1] != "/": args.data_directory += "/"
  resourceInfo['path'] = args.data_directory + "GRCh" + str(args.reference) + "/"

  # Check that the public-utils directory exists and add to the path so scripts from here can be used
  if args.utils_directory[-1] != "/": args.utils_directory += "/"
  if not os.path.exists(args.utils_directory): fail('public-utils directory could not be found')

  # Define the name of the resource description json file. Use the file provided on the command
  # line, and resort to the default file if not included. The version of the resource file is
  # unknown, so look for any json with the correct prefix and if there is more than one such
  # file, throw an error.
  if args.resource_json: resourceFilename = args.resource_json
  else:
    resourceFilename = args.data_directory + 'resources_GRCh' + str(args.reference) + '.json'
    resourceFiles    = glob.glob(resourceFilename)
    if len(resourceFiles) != 1: fail('There are zero, or more than one resource files for GRCh' + args.reference + ' in ' + args.data_directory)
    resourceFilename = resourceFiles[0]
  resourceInfo['json'] = resourceFilename

  # Check the family type is allowed
  if args.family_type not in allowedFamily: fail('The allowed family types (--family-type, -f) are: ' + ', '.join(allowedFamily))

  # If this is a singleton, no family based analyis should be performed
  resourceInfo['family_type'] = args.family_type
  resourceInfo['isFamily'] = False if args.family_type == 'singleton' else True

  # Check that the input files exist and have the correct extensions
  if not exists(args.input_vcf): fail('The vcf file could not be found')
  elif not args.input_vcf.endswith('.vcf.gz'): fail('The input vcf file (--vcf, -v) must be a compressed, indexed vcf and have the extension ".vcf.gz"')
  if resourceInfo['isFamily']:
    if not args.ped: fail('A ped file needs to specified (--ped, -e) for family type "' + str(args.family_type) + '"')
    if not exists(args.ped): fail('The ped file could not be found')
    elif not args.ped.endswith('.ped'): fail('The input ped file (--ped, -p) must have the extension ".ped"')

  # Return the info on resources
  return resourceInfo

# Build the path to allow importing additional modules
def buildPath(args):

  # Add public-utils and required directories to the path
  path.append(args.utils_directory)
  path.append(args.utils_directory + 'common_components')
  path.append(args.utils_directory + 'scripts')
  path.append(args.utils_directory + 'api_commands')

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
version = "1.0.0"
date    = str(date.today())

# The working directory where all created files are kept
workingDir = os.getcwd() + "/calypso_v" + version + "r"

## Store info on allowed values
allowedFamily     = ['singleton', 'duo', 'trio', 'quad', 'quintet']
allowedReferences = ['37', '38']

# Store information related to Mosaic
mosaicConfig = {}

if __name__ == "__main__":
  main()