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
import calypso_variant_filters as vfilt

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
  import api_sample_files as api_sf
  import api_sample_hpo_terms as api_shpo
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

  # Get the project reference
  reference = api_ps.getProjectReference(mosaicConfig, args.project_id) if not args.reference else args.reference
  if reference not in allowedReferences: fail('The specified reference (' + str(reference) + ') is not recognised. Allowed values are: ' + str(', '.join(allowedReferences)))
  print('  Using the reference:                  ', reference, sep = '')

  # Read the resources json file to identify all the resources that will be used in the annotation pipeline
  if args.resource_json: resourceInfo = res.checkResources(reference, args.data_directory, args.tools_directory, args.resource_json)
  else: resourceInfo = res.checkResources(reference, args.data_directory, args.tools_directory, args.data_directory + 'resources_' + str(reference) + '.json')
  resourceInfo = res.readResources(reference, rootPath, resourceInfo)

  # Define the tools to be used by Calypso
  resourceInfo = res.calypsoTools(resourceInfo)

  # Define paths to be used by Calypso
  setWorkingDir(resourceInfo['version'])

  # Determine the id of the proband, and the family structure
  args, proband, samples, familyType = sam.getProband(mosaicConfig, args, workingDir, api_s, api_ped)
  print('  Performing analysis on a family type: ', familyType, sep = '')
  if len(proband) == 0: print('  No proband could be determined:       ')
  elif len(proband) == 1: print('  Proband has the name:                 ', proband[0], sep = '')
  else:                 print('  Pedigree has multiple probands:       ', ', '.join(proband), sep = '')

  # Get the vcf file from Mosaic, if it can be determined
  if not args.input_vcf:
    vcfFiles = []
    for sample in proband:
      sampleVcfFiles = api_sf.getSampleVcfs(mosaicConfig, args.project_id, samples[sample]['mosaicId'])
      if len(sampleVcfFiles) == 1:
        uri = sampleVcfFiles[list(sampleVcfFiles.keys())[0]]['uri'][6:]
        if uri not in vcfFiles: vcfFiles.append(uri)
      else: fail('Mosaic has no, or multiple vcf files associated with it, so the file to use cannot be determined')
      if len(vcfFiles) != 1: fail('Mosaic has no, or multiple vcf files associated with it, so the file to use cannot be determined')
      args.input_vcf = vcfFiles[0]
  print('  Using the vcf file:                   ', args.input_vcf, sep = '')

  # Determine the order of the samples in the vcf header
  sam.getSampleOrder(resourceInfo, samples, args.input_vcf)

  # If the HPO terms are not specified on the command line, grab them from the project. If they are specified on the command
  # line, use these
  if not args.hpo:
    hpoTerms = []
    for sample in proband:
      sampleHpoTerms = api_shpo.getSampleHpoList(mosaicConfig, args.project_id, samples[sample]['mosaicId'])
      for hpo in sampleHpoTerms:
        if hpo not in hpoTerms: hpoTerms.append(hpo)
    args.hpo = ','.join(hpoTerms)
  print('  Using the HPO terms:                  ', args.hpo, sep = '')

  # Get all the project attributes in the target Mosaic project
  projectAttributes = mos.getProjectAttributes(mosaicConfig, args.project_id, api_pa)
  publicAttributes  = mos.getPublicProjectAttributes(mosaicConfig, api_pa)

  # Read the Mosaic json and validate its contents
  if not args.mosaic_json: args.mosaic_json = args.data_directory + 'resources_mosaic_' + reference + '.json'
  mosaicInfo         = mosr.readMosaicJson(args.mosaic_json, reference)
  projectAnnotations = mos.getProjectAnnotations(mosaicConfig, args.project_id, api_va)
  publicAnnotations  = mos.getPublicAnnotations(mosaicConfig, args.project_id, api_va)
  privateAnnotations, projectAnnotations = mos.createPrivateAnnotations(mosaicConfig, mosaicInfo['resources'], projectAnnotations, samples, args.project_id, api_va)
  mosr.checkPublicAnnotations(mosaicConfig, mosaicInfo['resources'], publicAnnotations, args.project_id, api_va)

  # Build the toml file for vcfanno. This defines all of the annotations that vcfanno needs to use
  # with the path to the required files
  tomlFilename = toml.buildToml(workingDir, resourceInfo)

  # Generate the bash script to run the annotation pipeline
  pipelineModifiers = bashScript.determinePipeline(familyType)
  bashFilename, bashFile = bashScript.openBashScript(workingDir)
  filteredVcf, clinvarVcf, rareVcf = bashScript.bashResources(resourceInfo, workingDir, bashFile, args.input_vcf, args.ped, tomlFilename, familyType, pipelineModifiers)
  bashScript.samplesFile(bashFile)
  bashScript.annotateVcf(resourceInfo, bashFile)
  #bashScript.annotateVcf(bashFile, familyType)#resourceInfo, pipelineModifiers)
  #bashScript.cleanedVcf(bashFile)
  #bashScript.annotateVcf(bashFile, resourceInfo, pipelineModifiers)
  bashScript.filterVariants(bashFile, samples, resourceInfo)
  #bashScript.rareDiseaseVariants(bashFile)
  bashScript.clinVarVariants(bashFile)

  # Perform filtering based on the family type. This filter will generate a small set of prioritized
  # variants based on the available family structure
  if familyType == 'trio': slivarTrio.annotateTrio(resourceInfo, bashFile)
  else: print('**** WARNING: Filtering for family type ' + str(familyType) + ' has not yet been implemented')
  bashScript.deleteFiles(args, pipelineModifiers, bashFile)

  # Process the filtered vcf file to extract the annotations to be uploaded to Mosaic
  print('# Generate tsv files to upload annotations to Mosaic', file = bashFile)
  uploadFilename, uploadFile = upload.openUploadAnnotationsFile(workingDir)

  # Define the output VCF files whose annotations need uploading to Mosaic
  annotateFiles = ['$FILTEREDVCF', '$CLINVARVCF']

  # Loop over all the resources to be uploaded to Mosaic
  for resource in mosaicInfo['resources']:
    if resource == 'hpo':

      # Only proceed with HPO terms if an HPO term string has been provided
      info     = mosaicInfo['resources']['hpo']
      fileInfo = resourceInfo['resources']['hpo']['file']
      if args.hpo: tsvFiles = anno.createHpoTsv(info, os.path.dirname(__file__) + '/components', args.config, args.hpo, args.project_id, fileInfo, args.utils_directory, args.tools_directory, bashFile, annotateFiles)
    else: tsvFiles = anno.createAnnotationTsv(mosaicInfo, resource, os.path.dirname(__file__) + '/components', reference, args.config, args.mosaic_json, args.tools_directory, bashFile, annotateFiles)
    for tsv in tsvFiles: upload.uploadAnnotations(args.utils_directory, tsv, mosaicInfo['resources'][resource]['project_id'], args.config, uploadFile)
  upload.closeUploadAnnotationsFile(uploadFilename, uploadFile)

  # Close the bash script
  bashScript.finishScript(bashFile, bashFilename, version)

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
  vfilt.readRequiredFilters(mosaicConfig, mosaicInfo, args.data_directory, samples, projectAnnotations, familyType, args.project_id, api_ps, api_vf)

  # Update Calypso attributes, e.g. the version that was run, the history of runs etc.
  mos.updateCalypsoAttributes(mosaicConfig, resourceInfo['version'], projectAttributes, publicAttributes, version, args.project_id, api_pa)

  # Generate scripts to upload filtered variants to Mosaic
  upload.uploadVariants(workingDir, args.utils_directory, args.config, args.project_id, filteredVcf, clinvarVcf, rareVcf)

  # Output a summary file listing the actions undertaken by Calypso with all version histories
  res.calypsoSummary(workingDir, version, resourceInfo, reference)
  print('Calypso pipeline version ', version, ' completed successfully', sep = '')

# Input options
def parseCommandLine():
  global version
  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  #parser.add_argument('--family_type', '-f', required = True, metavar = 'string', help = 'The familty structure. Allowed values: ' + ', '. join(allowedFamily))
  parser.add_argument('--data_directory', '-d', required = True, metavar = 'string', help = 'The path to the directory where the resources live')
  parser.add_argument('--tools_directory', '-s', required = False, metavar = 'string', help = 'The path to the directory where the tools to use live')
  parser.add_argument('--utils_directory', '-l', required = True, metavar = 'string', help = 'The path to the public-utils directory')
  parser.add_argument('--input_vcf', '-i', required = False, metavar = 'string', help = 'The input vcf file to annotate')
  parser.add_argument('--ped', '-e', required = False, metavar = 'string', help = 'The pedigree file for the family. Not required for singletons')
  parser.add_argument('--project_id', '-p', required = True, metavar = 'string', help = 'The project id that variants will be uploaded to')
  parser.add_argument('--config', '-c', required = True, metavar = "string", help = "The config file for Mosaic")

  # Optional pipeline arguments
  parser.add_argument('--reference', '-r', required = False, metavar = 'string', help = 'The reference genome to use. Allowed values: ' + ', '.join(allowedReferences))
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
