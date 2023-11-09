#!/usr/bin/python

from __future__ import print_function
from datetime import date
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
import tools_bcftools as bcftools

def main():
  global bcftoolsExe
  global slivarExe

  # Parse the command line
  args = parseCommandLine()

  # Define the executable bcftools command
  if args.tools_directory:
    if not args.tools_directory.endswith('/'): args.tools_directory = str(args.tools_directory) + '/'
  if args.tools_directory: bcftoolsExe = args.tools_directory + 'bcftools/bcftools'
  else: bcftoolsExe = 'bcftools'

  # ...and the slivar executable command
  if args.tools_directory: slivarExe = args.tools_directory + 'slivar'
  else: slivarExe = 'slivar'

  # Read the json describing the variants to report
  newSettings, removedSettings, updatedSettings = readSettings(args.json_file)

  # Check that the vcf files (the previous and the current) exist and extract the dates of these
  # files
  previousDate, currentDate = checkVcfs(args.previous_vcf, args.current_vcf)
  output = 'diff_' + str(previousDate) + '_' + str(currentDate) + '.vcf.gz'

  # Generete the intsersection and complements of the two clinVar vcf files
  command = bcftools.isec(bcftoolsExe, args.previous_vcf, args.current_vcf, str(workingDir) + 'diff').split(' ')
  process = subprocess.Popen(command, stdout = subprocess.PIPE)
  process.wait()

  # Process the vcf files containing the new clinVar variants (0000.vcf.gz) and the removed (0001.vcf.gz).
  # The output vcf will only contain variants specified in newSettings as variants to report
  processUniqueVariants(newSettings, 'diff/0000.vcf.gz', 'new', str(workingDir) + 'diff/new.vcf')
  processUniqueVariants(removedSettings, 'diff/0001.vcf.gz', 'removed', str(workingDir) + 'diff/removed.vcf')

  # To identify updated clinVar variants, use the 0002.vcf.gz and 0003.vcf.gz files. These files contain the
  # variants shared by the two clinVar files, but with the annotations from the previous (0002) and the current
  # (0003) files.
  processupdatedVariants(updatedSettings, args.convert_annotations, 'updated', str(workingDir) + 'diff/updated.vcf')

  # Merge the output files containing the new, removed, and updated variants to report
  mergeFiles('diff/new.vcf', 'diff/removed.vcf', 'diff/updated.vcf', output)

# Merge the files containing the variants to report, sort, compress and index
def mergeFiles(new, removed, updated, output):
  global workingDir

  #

# Input options
def parseCommandLine():
  global version

  parser = argparse.ArgumentParser(description='Process the command line')

  # Mosaic arguments
  parser.add_argument('--config', '-c', required = False, metavar = "string", help = "The config file for Mosaic")
  parser.add_argument('--token', '-t', required = False, metavar = "string", help = "The Mosaic authorization token")
  parser.add_argument('--url', '-u', required = False, metavar = "string", help = "The base url for Mosaic")

  # The project id
  #parser.add_argument('--project_id', '-p', required = True, metavar = "string", help = "The project id that variants will be uploaded to")

  # Directories where required tools live
  #parser.add_argument('--utils_directory', '-l', required = True, metavar = 'string', help = 'The path to the public-utils directory')
  parser.add_argument('--tools_directory', '-s', required = False, metavar = 'string', help = 'The path to the tools directory')

  # A json file describing the variants to report is required
  parser.add_argument('--json_file', '-j', required = True, metavar = 'string', help = 'The json file describing which variants to report')

  # The previous and current clinVar vcf files to compare
  parser.add_argument('--previous_vcf', '-r', required = True, metavar = 'string', help = 'The previous clinVar vcf file')
  parser.add_argument('--current_vcf', '-e', required = True, metavar = 'string', help = 'The current clinVar vcf file')

  # The file that defines how to rename annotations
  parser.add_argument('--convert_annotations', '-o', required = True, metavar = 'string', help = 'The file with description of how to convert annotations')

  # Version
  parser.add_argument('--version', '-v', action="version", version='HPO prioritization pipeline version: ' + str(version))

  return parser.parse_args()

# Read the json file describing which new, removed, or updated clinVar variants to report
def readSettings(settingsJson):

  # Open the json file
  try: settingsFile = open(settingsJson, 'r')
  except: fail('Unable to open json settings file: ' + str(settingsJson))

  # Load this as a json file
  try: settings = json.load(settingsFile)
  except: fail('Unable to read settings json file. Ensure this is a valid json file.')

  # Ensure that the "new", "removed", and "updated" sections are included to define how to handle these three sections
  if 'new' not in settings: fail('The section "new" is not included in the settings json file')
  if 'removed' not in settings: fail('The section "removed" is not included in the settings json file')
  if 'updated' not in settings: fail('The section "updated" is not included in the settings json file')

  # Store the information for each
  newSettings     = settings['new']
  removedSettings = settings['removed']
  updatedSettings = settings['updated']

  # Return the settings
  return newSettings, removedSettings, updatedSettings

# Check that the vcfs to compare exist
def checkVcfs(previous, current):

  # Check that the files exist
  if not exists(previous): fail('Could not open file: ' + str(previous))
  if not exists(current): fail('Could not open file: ' + str(current))

  # Get and return the dates of the clinVar files
  previousDate = previous.split('/')[-1].split('.')[0].split('_')[-1]
  currentDate  = current.split('/')[-1].split('.')[0].split('_')[-1]
  return previousDate, currentDate

# Process the vcf files with unique variants
def processUniqueVariants(settings, vcf, fileId, output):
  global bcftoolsExe

  # Create the output file using the header of the input file
  command = bcftools.createHeaderVcf(bcftoolsExe, vcf, output).split(' ')
  process = subprocess.Popen(command, stdout = subprocess.PIPE)
  process.wait()

  # Open the output file
  outputFile = open(output, 'w')

  # Loop over all lines in the vcf file having extracted the CLNSIG annotation
  for line in os.popen(bcftools.query(bcftoolsExe, vcf, ['CLNSIG'])).readlines():
    value = line.rstrip().split('\t')[5]
    if value != '.': 

      # If the clinical significance has multiple values (separated by a ',' or a '|'), break the value
      # into its constituent pieces
      if ',' in value: significance = value.split(',')
      elif '|' in value: significance = value.split('|')
      else: significance = [value]

      # Loop over the available values. If any of them are listed as a significance to report, write out the line
      reportVariant = False
      for significanceValue in significance:
        try: 
          if settings[significanceValue]:
            reportVariant = True
            break
        except:
          details = line.rstrip().split()
          fail('Unknown significance (' + str(significanceValue) + ') in ' + str(fileId) + ' file at position: ' + str(details[0]) + ':' + str(details[1]))
      if reportVariant: print(line.rstrip(), file = outputFile)

  # Close the output file
  outputFile.close()

# Process the variants shared by both ClinVar files
def processUpdatedVariants(settings, convertAnnotations, fileId, output):
  global bcftoolsExe
  global slivarExe
  global mendelianSignificance
  global workingDir

  # Check that the file describing how to convert annotations exists
  if not exists(convertAnnotations): fail('Could not open file ' + str(convertAnnotations))

  # Open the output file
  outputFile = open(output, 'w')

  # Extract the CLNSIG annotations from the previous ClinVar vcf and rename CLNSIG to PCLNSIG
  command = bcftools.extractAnnotations(bcftoolsExe, ['CLNSIG'], convertAnnotations, str(workingDir) + 'diff/0002.vcf.gz', str(workingDir) + 'diff/previous_shared.vcf.gz').split(' ')
  process = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  stdout, stderr = process.communicate()
  process.wait()

  # Extract the CLNSIG annotations from the current ClinVar vcf
  command = bcftools.extractAnnotations(bcftoolsExe, ['CLNSIG'], False, str(workingDir) + 'diff/0003.vcf.gz', str(workingDir) + 'diff/current_shared.vcf.gz').split(' ')
  process = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  stdout, stderr = process.communicate()
  process.wait()

  # Generate a new vcf with the annotations from both the previous and current ClinVar vcf files
  command = bcftools.annotateWithVcf(bcftoolsExe, ['PCLNSIG'], str(workingDir) + 'diff/previous_shared.vcf.gz', str(workingDir) + 'diff/current_shared.vcf.gz', False)
  command += ' | ' + slivarExe + ' expr --vcf - --info \'("PCLNSIG" in INFO && "CLNSIG" in INFO && INFO.PCLNSIG != INFO.CLNSIG)\' --pass-only'
  command += ' | ' + bcftools.compressAndIndex(bcftoolsExe, '-', str(workingDir) + 'diff/shared.vcf.gz', 'z')
  process = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  stdout, stderr = process.communicate()
  process.wait()

  # Now loop over the shared.vcf.gz file and determine which variants to report
  for line in os.popen(bcftools.query(bcftoolsExe, str(workingDir) + 'diff/shared.vcf.gz', ['PCLNSIG', 'CLNSIG'])).readlines():
    fields = line.rstrip().split('\t')

    # If the clinical significance has multiple values (separated by a ',' or a '|'), break the value
    # into its constituent pieces
    annotations = fields[5]
    if ',' in annotations: previousAnnotations = annotations.split(',')
    elif '|' in annotations: previousAnnotations = annotations.split('|')
    else: previousAnnotations = [annotations]

    annotations = fields[6]
    if ',' in annotations: currentAnnotations = annotations.split(',')
    elif '|' in annotations: currentAnnotations = annotations.split('|')
    else: currentAnnotations = [annotations]

    # Get the Mendelian signicifance value. This is one of the terms Pathogenic, Likely_pathogenic, Uncertain_significance,
    # Likely_benign, or Benign. Also Conflicting interpretations can be accepted. There can only be one of these values
    previousSignificance = 'Null'
    for annotation in previousAnnotations:
      if annotation in mendelianSignificance:
        previousSignificance = annotation
        break
    currentSignificance = 'Null'
    for annotation in currentAnnotations:
      if annotation in mendelianSignificance:
        currentSignificance = annotation
        break

    # Check if the listed change should be reported
    reportVariant = False
    if previousSignificance not in settings: 
      fail('Unknown significance (' + str(previousSignificance) + ') in ' + str(fileId) + ' file at position: ' + str(fields[0]) + ':' + str(fields[1]))
    if currentSignificance not in settings[previousSignificance]:
      fail('Unknown significance (' + str(currentSignificance) + ') in the ' + str(previousSignificance) + ' settings for the ' + str(fileId) + ' file at position: ' + str(fields[0]) + ':' + str(fields[1]))
    if settings[previousSignificance][currentSignificance]: reportVariant = True

    # Handle any of the additional significance annotations, e.g. low_penetrance recommandations from ClinGen, or other
    # allowable terms (e.g. drug_response)
##########
##########
########## COMPLETE
##########
##########

    # If the variant is to be reported, write it out
    if reportVariant: print(fields[0], fields[1], fields[2], fields[3], fields[4], 'CLNSIG=' + str(fields[6]), sep = '\t', file = outputFile)

  # Close the output file
  outputFile.close()

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = '')
  exit(1)

# Initialise global variables

# Get the directory where the files will be generated
workingDir = os.getcwd() + '/'

# Version
version = "0.0.1"
date    = str(date.today())

# The bcftools executable command
bcftoolsExe = False

# The allowed CLNSIG values for Mendelian variants. The term Null is allowed here. This is to handle cases
# where the CLNSIG is e.g. drug_response. There is still a CLNSIG value, but it is not one of the ACMG/AMP
# terms for a Mendelian variant. In the case, the significance can be listed as Null, so a variant that
# changes from Pathogenic to drug_reponse would appear as Pathogenic > Null and could be reported
mendelianSignificance = []
mendelianSignificance.append('Pathogenic')
mendelianSignificance.append('Pathogenic/Likely_pathogenic')
mendelianSignificance.append('Likely_pathogenic')
mendelianSignificance.append('Uncertain_significance')
mendelianSignificance.append('Conflicting_interpretations_of_pathogenicity')
mendelianSignificance.append('Likely_benign')
mendelianSignificance.append('Benign/Likely_benign')
mendelianSignificance.append('Benign')
mendelianSignificance.append('Null')

if __name__ == "__main__":
  main()
