from os.path import exists

import argparse
import os
import subprocess
import sys

# Add the path of the common functions and import them
from sys import path
path.append(os.path.dirname(__file__) + '/../components')
import calypso_path as cpath
import tools_bcftools as bcftools

def main():
  global bcftoolsExe
  global slivarExe

  # Parse the command line
  args = parseCommandLine()

  # Check the supplied parameters are as expected, then expand the path to enable scripts and api commands
  # to be accessed for Calypso
  rootPath = os.path.dirname(__file__)
  if not args.tools_directory.endswith('/'): args.tools_directory = str(args.tools_directory) + '/'
  sys.path = cpath.buildPath(sys.path, args.utils_directory)

  # Import additional tools
  import mosaic_config
  import api_samples as api_s
  import api_sample_attributes as api_sa
  import api_sample_hpo_terms as api_sh

  # Parse the Mosaic config file to get the token and url for the api calls
  mosaicRequired = {'MOSAIC_TOKEN': {'value': args.token, 'desc': 'An access token', 'long': '--token', 'short': '-t'},
                    'MOSAIC_URL': {'value': args.url, 'desc': 'The api url', 'long': '--url', 'short': '-u'},
                    'MOSAIC_ATTRIBUTES_PROJECT_ID': {'value': args.attributes_project, 'desc': 'The public attributes project id', 'long': '--attributes_project', 'short': '-a'}}
  mosaicConfig   = mosaic_config.mosaicConfigFile(args.config)
  mosaicConfig   = mosaic_config.commandLineArguments(mosaicConfig, mosaicRequired)

  # Define the executable bcftools command
  bcftoolsExe = args.tools_directory + 'bcftools/bcftools'

  # ...and the slivar executable command
  slivarExe = args.tools_directory + 'slivar'

  # Get the Mosaic sample id for the provided sample
  sampleAttributes = api_sa.getSampleAttributesWValues(mosaicConfig, args.project_id)
  isProband = False
  probandId = False
  for attribute in sampleAttributes:
    if attribute['name'] == 'Relation':
      for sampleInfo in attribute['values']:
        if sampleInfo['value'] == 'Proband':
          if isProband: fail('Multiple probands in project with id ' + str(args.project_id))
          isProband = True
          probandId = sampleInfo['sample_id']
          break
  if not isProband: fail('No proband found for project with id ' + str(args.project_id))
  sampleIds   = api_s.getSamplesDictIdName(mosaicConfig, args.project_id)
  probandName = False
  for sampleId in sampleIds:
    if sampleId == probandId: probandName = sampleIds[sampleId]
    break

  # Get required case information from Mosaic
  hpoTerms = api_sh.getSampleHpo(mosaicConfig, args.project_id, sampleId)
  print(hpoTerms)

# Input options
def parseCommandLine():
  global version

  parser = argparse.ArgumentParser(description='Process the command line')

  # Mosaic arguments
  parser.add_argument('--config', '-c', required = False, metavar = "string", help = "The config file for Mosaic")
  parser.add_argument('--token', '-t', required = False, metavar = "string", help = "The Mosaic authorization token")
  parser.add_argument('--url', '-u', required = False, metavar = "string", help = "The base url for Mosaic")
  parser.add_argument('--attributes_project', '-a', required = False, metavar = "integer", help = "The Mosaic project id that contains public attributes")

  # The project id
  parser.add_argument('--project_id', '-p', required = True, metavar = "string", help = "The project id that variants will be uploaded to")

  # Directories where required tools live
  parser.add_argument('--utils_directory', '-l', required = True, metavar = 'string', help = 'The path to the public-utils directory')
  parser.add_argument('--tools_directory', '-s', required = True, metavar = 'string', help = 'The path to the tools directory')

  # A json file describing the variants to report is required
  #parser.add_argument('--json_file', '-j', required = True, metavar = 'string', help = 'The json file describing which variants to report')

  # The previous and current clinVar vcf files to compare
  parser.add_argument('--vcf', '-i', required = True, metavar = 'string', help = 'The vcf file containing the differences between two versions')

  # Version
  parser.add_argument('--version', '-v', action="version", version='HPO prioritization pipeline version: ' + str(version))

  return parser.parse_args()

# Check that the vcfs to compare exist
def checkVcfs(previous, current):

  # Check that the files exist
  if not exists(previous): fail('Could not open file: ' + str(previous))
  if not exists(current): fail('Could not open file: ' + str(current))

  # Get and return the dates of the clinVar files
  previousDate = previous.split('/')[-1].split('.')[0].split('_')[-1]
  currentDate  = current.split('/')[-1].split('.')[0].split('_')[-1]
  return previousDate, currentDate

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = '')
  exit(1)

# Initialise global variables

# Get the directory where the files will be generated
workingDir = os.getcwd() + '/'

# Version
version = "0.0.1"

# The bcftools executable command
bcftoolsExe = False

if __name__ == "__main__":
  main()
