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

  # Import the api client
  path.append(args.api_client)
  from mosaic import Mosaic, Project, Store
  store   = Store(config_file = args.config)
  mosaic  = Mosaic(config_file = args.config)

  # Create a project object for the defined project
  project = mosaic.get_project(args.project_id)

  # Check the supplied parameters are as expected, then expand the path to enable scripts and api commands
  # to be accessed for Calypso
  rootPath = os.path.dirname(__file__)
  if not args.tools_directory.endswith('/'): args.tools_directory = str(args.tools_directory) + '/'

  # Define the executable bcftools command
  bcftoolsExe = args.tools_directory + 'bcftools/bcftools'

  # ...and the slivar executable command
  slivarExe = args.tools_directory + 'slivar'

  # Get the Mosaic sample id for the provided sample
  isProband = False
  probandId = False
  samples   = project.get_samples()
  for sample in samples:
    sampleId = sample['id']
    data     = project.get_attributes_for_sample(sampleId)
    for attribute in data:
      if attribute['name'] == 'Relation':

        # Loop over the values and see if this is the proband
        for sampleValues in attribute['values']:
          if sampleValues['value'] == 'Proband':
            if isProband: fail('Multiple probands in project with id ' + str(args.project_id))
            isProband   = True
            probandId   = sampleValues['sample_id']
            probandName = sample['name']
            break

  # If no proband is found, fail
  if not isProband: fail('No proband found for project with id ' + str(args.project_id))

  # Get the HPO terms for the proband
  hpoTerms = {}
  for hpoTermInfo in project.get_sample_hpo_terms(probandId):
    hpoTerms[hpoTermInfo['hpo_id']] = hpoTermInfo['label']

  # Get a list of genes associated with the HPO terms using the HPO phenotype_to_gene info
  hpoFile = '/scratch/ucgd/lustre-work/marth/marth-projects/calypso/reannotation/data/GRCh38/hpo/phenotype_to_genes.txt'
  hpoGenes = open(hpoFile, 'r')
  genes    = []
  for line in hpoGenes.readlines():
    fields = line.rstrip().split('\t')
    term   = fields[0]
    gene   = fields[3]
    if term in hpoTerms:
      if gene not in genes: genes.append(gene)
  hpoGenes.close()
  print(genes)
  print(len(hpoTerm))
  print(len(genes))

# Input options
def parseCommandLine():
  global version

  parser = argparse.ArgumentParser(description='Process the command line')

  # Mosaic arguments
  parser.add_argument('--config', '-c', required = True, metavar = 'string', help = 'The config file for Mosaic')
  parser.add_argument('--api_client', '-a', required = True, metavar = 'string', help = 'The directory where the Python api wrapper lives')

  # The project id
  parser.add_argument('--project_id', '-p', required = True, metavar = "string", help = "The project id that variants will be uploaded to")

  # Directories where required tools live
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
