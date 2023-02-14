#!/usr/bin/python

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
  checkArguments(args)
  sys.path = cpath.buildPath(sys.path, args.utils_directory)

  # Read the resources json file to identify all the resources that will be used in the annotation pipeline
  if args.resource_json: resourceInfo = res.checkResources(args.reference, args.data_directory, args.resource_json)
  else: resourceInfo = res.checkResources(args.reference, args.data_directory, args.data_directory + 'resources_GRCh' + str(args.reference) + '.json')
  resourceInfo = res.readResources(args.reference, resourceInfo)
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
  import add_variant_filters as pu_avf

  # Parse the Mosaic config file to get the token and url for the api calls
  mosaicRequired = {'MOSAIC_TOKEN': {'value': args.token, 'desc': 'An access token', 'long': '--token', 'short': '-t'},
                    'MOSAIC_URL': {'value': args.url, 'desc': 'The api url', 'long': '--url', 'short': '-u'},
                    'MOSAIC_ATTRIBUTES_PROJECT_ID': {'value': args.attributes_project, 'desc': 'The public attribtes project id', 'long': '--attributes_project', 'short': '-a'}}
  mosaicConfig   = mosaic_config.mosaicConfigFile(args.config)
  mosaicConfig   = mosaic_config.commandLineArguments(mosaicConfig, mosaicRequired)

  # Determine the id of the proband, and the order of the samples in the vcf header
  proband, samples = sam.getProband(mosaicConfig, args.ped, args.family_type, args.project_id, api_s)
  sam.getSampleOrder(samples, args.input_vcf)

  # Get all the project attributes in the target Mosaic project
  projectAttributes = mos.getProjectAttributes(mosaicConfig, args.project_id, api_pa)
  publicAttributes  = mos.getPublicProjectAttributes(mosaicConfig, api_pa)

  # Get the value for the Calypso resource version attribute. If Calypso has been run previously, this will indicate the last resources
  # that were used
  #updatedResources        = []
  #previousResourceVersion = mos.getPreviousResourceVersion(projectAttributes)
  #if previousResourceVersion and str(previousResourceVersion) != str(resourceInfo['version']):
  #  previousResourceInfo = res.checkResources(args.reference, args.data_directory, args.data_directory + 'resources_archive/resources_GRCh' + str(args.reference) + '_v' + str(previousResourceVersion) + '.json')
  #  previousResourceInfo = res.readResources(args.reference, previousResourceInfo)
  #  updatedResources     = rean.compare(previousResourceInfo, resourceInfo)
#
#    # If Calypso was updated, see if the difference files for the original and updated ClinVar resources exist. If not, create them
#    if 'ClinVar' in updatedResources:
#      cvDir   = str(args.data_directory) + 'GRCh' + str(args.reference) + '/clinvar/'
#      cvDates = str(updatedResources['ClinVar']['previousVersion']) + '_' + str(updatedResources['ClinVar']['currentVersion'])
#      compDir = str(cvDir) + str(cvDates)
#      if not exists(compDir): cv.compare(compDir, previousResourceInfo['resources']['ClinVar']['file'], resourceInfo['resources']['ClinVar']['file'])

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
  filteredVcf, clinvarVcf, rareVcf = bashScript.bashResources(resourceInfo, workingDir, bashFile, args.input_vcf, args.ped, tomlFilename)
  bashScript.samplesFile(bashFile)
  bashScript.cleanedVcf(bashFile)
  bashScript.annotateVcf(bashFile, resourceInfo, args.family_type)
  bashScript.filterVariants(bashFile, proband, resourceInfo)
  bashScript.rareDiseaseVariants(bashFile)
  bashScript.clinVarVariants(bashFile)
  bashScript.deleteFiles(bashFile)

  # Process the filtered vcf file to extract the annotations to be uploaded to Mosaic
  print('# Generate tsv files to upload annotations to Mosaic', file = bashFile)
  uploadFilename, uploadFile = upload.openUploadAnnotationsFile(workingDir)

  # Define the output VCF files whose annotations need uploading to Mosaic
  annotateFiles = ['$FILTEREDVCF', '$CLINVARVCF']

  # Loop over all the resources to be uploaded to Mosaic
  for resource in mosaicInfo['resources']:
    if resource == 'hpo':

      # Only proceed with HPO terms if an HPO term string has been provided
      if args.hpo: tsvFiles = anno.createHpoTsv(mosaicInfo['resources']['hpo'], os.path.dirname(__file__) + '/components', args.config, args.hpo, args.project_id, resourceInfo['resources']['hpo']['file'], args.utils_directory, bashFile, annotateFiles)

    else: tsvFiles = anno.createAnnotationTsv(mosaicInfo, resource, os.path.dirname(__file__) + '/components', args.reference, args.config, args.mosaic_json, bashFile, annotateFiles)
    for tsv in tsvFiles: upload.uploadAnnotations(args.utils_directory, tsv, mosaicInfo['resources'][resource]['project_id'], args.config, uploadFile)
  upload.closeUploadAnnotationsFile(uploadFilename, uploadFile)

  # If this is a reannotation of the sample and there has been a change to ClinVar, perform a check for ClinVar significance
  # changes based on the difference files.
  #if 'ClinVar' in updatedResources: rean.clinvarUpdates(cvDates, compDir, bashFile)

  # Close the bash script
  bashScript.finishScript(bashFile, bashFilename, version)

  # Prepare the target project. This includes:
  # 1. Remove any unnecessary annotations from the project
  # 2. Import public annotations into the project
  # 3. Set the default annotations for the project
  # 4. Set up required variant filters
  # 5. Add and set project attributes
  if 'remove' in mosaicInfo: mos.removeAnnotations(mosaicConfig, mosaicInfo['remove'], projectAnnotations, args.project_id, api_va)
  mos.importAnnotations(mosaicConfig, mosaicInfo['resources'], projectAnnotations, publicAnnotations, args.project_id, api_va)
  mos.defaultAnnotations(mosaicConfig, mosaicInfo['defaultAnnotations'], publicAnnotations, privateAnnotations, args.project_id, api_ps)

############
############
############ SORT OUT NEW VARIANT FILTERS
############
############
  #pu_avf.addVariantFilters(mosaicConfig, args.filter_json, args.project_id, sampleIds, uids)
  #vfilt.variantFilters(mosaicConfig, mosaicInfo, rootPath, samples, privateAnnotations, args.project_id, api_vf)
  mos.updateCalypsoAttributes(mosaicConfig, resourceInfo['version'], projectAttributes, publicAttributes, version, args.project_id, api_pa)

  # Generate scripts to upload filtered variants to Mosaic
  upload.uploadVariants(workingDir, args.utils_directory, args.config, args.project_id, filteredVcf, clinvarVcf, rareVcf)

  # Output a summary file listing the actions undertaken by Calypso with all version histories
  res.calypsoSummary(workingDir, version, resourceInfo, args.reference)
  print('Calypso pipeline version ', version, ' completed successfully', sep = '')

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
  global allowedReferences
  global allowedFamily

  # Check the reference is allowed
  if args.reference not in allowedReferences:
    message = 'The allowed references (--reference, -r) are:'
    for ref in allowedReferences: message += '\n  ' + str(ref)
    fail(message)

  # Ensure the data path ends with a "/", then add the reference directory
  if args.data_directory[-1] != "/": args.data_directory += "/"

  # Check that the public-utils directory exists and add to the path so scripts from here can be used
  if args.utils_directory[-1] != "/": args.utils_directory += "/"
  if not os.path.exists(args.utils_directory): fail('public-utils directory could not be found')

  # Check the family type is allowed
  if args.family_type not in allowedFamily: fail('The allowed family types (--family-type, -f) are: ' + ', '.join(allowedFamily))

  # Check that the input files exist and have the correct extensions
  if not exists(args.input_vcf): fail('The vcf file could not be found')
  elif not args.input_vcf.endswith('.vcf.gz'): fail('The input vcf file (--vcf, -v) must be a compressed, indexed vcf and have the extension ".vcf.gz"')
  if args.family_type != 'singleton':
    if not args.ped: fail('A ped file needs to specified (--ped, -e) for family type "' + str(args.family_type) + '"')
    if not exists(args.ped): fail('The ped file could not be found')
    elif not args.ped.endswith('.ped'): fail('The input ped file (--ped, -p) must have the extension ".ped"')

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
version = "1.1.0"
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
