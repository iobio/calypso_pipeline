from os.path import exists

import argparse
import os
import json
import shutil
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
  if not args.tools_directory.endswith('/'): args.tools_directory = str(args.tools_directory) + '/'
  bcftoolsExe = args.tools_directory + 'bcftools/bcftools'

  # ...and the slivar executable command
  slivarExe = args.tools_directory + 'slivar'

  # Read the json describing the variants to report
  newSettings, removedSettings, updatedSettings = readSettings(args.json_file)

  # Check that the vcf files (the previous and the current) exist and extract the dates of these
  # files
  previousDate, currentDate = checkVcfs(args.previous_vcf, args.current_vcf)
  output = 'clinvar_differences_' + str(previousDate) + '_' + str(currentDate) + '.vcf.gz'

  # Generate the intsersection and complements of the two clinVar vcf files
  command = bcftools.isec(bcftoolsExe, args.previous_vcf, args.current_vcf, str(workingDir) + 'diff').split(' ')
  process = subprocess.Popen(command, stdout = subprocess.PIPE)
  process.wait()

  # Process the vcf files containing the new clinVar variants (0000.vcf.gz) and the removed (0001.vcf.gz).
  # The output vcf will only contain variants specified in newSettings as variants to report
  processUniqueVariants(newSettings, workingDir, '0000.vcf.gz', 'new.vcf.gz', 'new', ['NEW'])
  processUniqueVariants(removedSettings, workingDir, '0001.vcf.gz', 'removed.vcf.gz', 'removed', ['REM'])

  # To identify updated clinVar variants, use the 0002.vcf.gz and 0003.vcf.gz files. These files contain the
  # variants shared by the two clinVar files, but with the annotations from the previous (0002) and the current
  # (0003) files.
  processUpdatedVariants(updatedSettings, 'updated', str(workingDir) + 'diff/updated.vcf')

  # Merge the output files containing the new, removed, and updated variants to report
  mergeFiles('diff/new.vcf.gz', 'diff/removed.vcf.gz', 'diff/updated.vcf.gz', output)

  # Remove the 'diff' directory with all the intermediate files
  success = shutil.rmtree(str(workingDir) + 'diff')

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
  parser.add_argument('--tools_directory', '-s', required = True, metavar = 'string', help = 'The path to the tools directory')

  # A json file describing the variants to report is required
  parser.add_argument('--json_file', '-j', required = True, metavar = 'string', help = 'The json file describing which variants to report')

  # The previous and current clinVar vcf files to compare
  parser.add_argument('--previous_vcf', '-r', required = True, metavar = 'string', help = 'The previous clinVar vcf file')
  parser.add_argument('--current_vcf', '-e', required = True, metavar = 'string', help = 'The current clinVar vcf file')

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
def processUniqueVariants(settings, workingDir, inputVcf, outputVcf, fileId, addedAnnotations):
  global bcftoolsExe
  global stars

  # Add the path to the vcf files
  inputVcf  = str(workingDir) + 'diff/' + str(inputVcf)
  outputVcf = str(workingDir) + 'diff/' + str(outputVcf)

  # Write new header lines to a file
  headerName = str(workingDir) + 'diff/header.vcf'
  header     = open(headerName, 'w')
  for annotation in addedAnnotations: print('##INFO=<ID=' + str(annotation) + ',Number=1,Type=String,Description=Tag describing update in ClinVar>', sep = '', file = header)
  header.close()
  
  # Run the bcftools annotate command to extract the required annotations
  command = bcftools.extractAnnotations(bcftoolsExe, ['CLNSIG', 'CLNREVSTAT'], addedAnnotations, headerName, False, inputVcf, outputVcf)
  process = subprocess.Popen(command.split(' '), stdout = subprocess.PIPE)
  process.wait()

#  # Create the output file using the header of the input file
#  command = bcftools.createHeaderVcf(bcftoolsExe, vcf, output).split(' ')
#  process = subprocess.Popen(command, stdout = subprocess.PIPE)
#  process.wait()
#
#  # Open the output file
#  outputFile = open(output, 'a')
#
#  # Loop over all lines in the vcf file having extracted the CLNSIG annotation
#  for line in os.popen(bcftools.query(bcftoolsExe, vcf, ['CLNSIG', 'CLNREVSTAT'])).readlines():
#    fields = line.rstrip().split('\t')
#    clnsig = fields[5]
#
#    # Replace the CLNREVSTAT text with the star value
#    try: fields[6] = stars[fields[6]]
#    except: fail('CLNREVSTAT has unknown value: ' + str(fields[6]))
#    if clnsig != '.': 
#
#      # If the clinical significance has multiple values (separated by a ',' or a '|'), break the value
#      # into its constituent pieces
#      if ',' in clnsig: significance = clnsig.split(',')
#      elif '|' in clnsig: significance = clnsig.split('|')
#      else: significance = [clnsig]
#
#      # Loop over the available values. If any of them are listed as a significance to report, write out the line
#      reportVariant = False
#      for significanceValue in significance:
#
#        # The settings contains all the terms that should be reported. If any of these terms are present in the
#        # significance, the variant should be reported. For example, Likely_pathogenic could be a term in the
#        # settings which would allow Likely_pathogenic/Likely_risk_allele to be reported.
#        for reportableSig in settings:
#          if str(reportableSig.lower()) in str(significanceValue.lower()):
#            reportVariant = True
#            break
#      if reportVariant: print('\t'.join(fields), file = outputFile)
#
#  # Close the output file
#  outputFile.close()

# Process the variants shared by both ClinVar files
def processUpdatedVariants(settings, fileId, output):
  global auditing
  global bcftoolsExe
  global slivarExe
  global ignoredPclnsig
  global stars
  global workingDir

  # Create a file to use to convert the CLNSIG and CLNREVSTAT annotations to PCLNSIG and PCLNREVSTAT for the previous
  # ClinVar file
  convName = str(workingDir) + 'diff/convert_annotations.txt'
  convFile = open(convName, 'w')
  print('INFO/CLNSIG\tPCLNSIG', file = convFile)
  print('INFO/CLNREVSTAT\tPCLNREVSTAT', file = convFile)
  convFile.close()

  # Extract the CLNSIG annotations from the previous ClinVar vcf and rename CLNSIG to PCLNSIG
  command = bcftools.extractAnnotations(bcftoolsExe, ['CLNSIG', 'CLNREVSTAT'], False, False, convName, str(workingDir) + 'diff/0002.vcf.gz', str(workingDir) + 'diff/previous_shared.vcf.gz').split(' ')
  process = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  stdout, stderr = process.communicate()
  process.wait()

  # Extract the CLNSIG annotations from the current ClinVar vcf
  command = bcftools.extractAnnotations(bcftoolsExe, ['CLNSIG', 'CLNREVSTAT'], False, False, False, str(workingDir) + 'diff/0003.vcf.gz', str(workingDir) + 'diff/current_shared.vcf.gz').split(' ')
  process = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  stdout, stderr = process.communicate()
  process.wait()

  # Generate a new vcf with the annotations from both the previous and current ClinVar vcf files
  command = bcftools.annotateWithVcf(bcftoolsExe, ['PCLNSIG', 'PCLNREVSTAT'], str(workingDir) + 'diff/previous_shared.vcf.gz', str(workingDir) + 'diff/current_shared.vcf.gz', False)
  command += ' | ' + slivarExe + ' expr --vcf - --info \'("PCLNSIG" in INFO && "CLNSIG" in INFO && INFO.PCLNSIG != INFO.CLNSIG)\' --pass-only'
  command += ' | ' + bcftools.compressAndIndex(bcftoolsExe, '-', str(workingDir) + 'diff/shared.vcf.gz', 'z')
  process = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  stdout, stderr = process.communicate()
  process.wait()

  # Create the output file using the header of the input file
  command = bcftools.createHeaderVcf(bcftoolsExe, str(workingDir) + 'diff/shared.vcf.gz', output)
  process = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  process.wait()

  # Open the output file
  outputFile = open(output, 'a')

  # Now loop over the shared.vcf.gz file and determine which variants to report
  totalVariants = 0
  for line in os.popen(bcftools.queryVcf(bcftoolsExe, str(workingDir) + 'diff/shared.vcf.gz', ['PCLNSIG', 'CLNSIG', 'PCLNREVSTAT', 'CLNREVSTAT'])).readlines():
    fields        = line.rstrip().split('\t')
    reportVariant = False
    totalVariants += 1

    # Convert the CLNREVSTAT fields to the number of stars
    try: pclnrev = stars[fields[9]]
    except: fail('Unknown CLNREVSTAT value: ' + str(fields[9]))
    try: clnrev = stars[fields[10]]
    except: fail('Unknown CLNREVSTAT value: ' + str(fields[10]))

    # If the clinical significance has multiple values (separated by a ',', a '/' or a '|'), break the value
    # into its constituent pieces
    if ',' in fields[7]: pclnsigs = fields[7].split(',')
    elif '/' in fields[7]: pclnsigs = fields[7].split('/')
    elif '|' in fields[7]: pclnsigs = fields[7].split('|')
    else: pclnsigs = [fields[7]]

    if ',' in fields[8]: clnsigs = fields[8].split(',')
    elif '/' in fields[8]: clnsigs = fields[8].split('/')
    elif '|' in fields[8]: clnsigs = fields[8].split('|')
    else: clnsigs = [fields[8]]

    # Loop over the values in the previous significance
    for pclnsig in pclnsigs:
      if pclnsig in settings:

        # If 'not' occurs in the settings, this means that the variant should be reported if the CLNSIG
        # in the new ClinVar vcf is 'not' any of these values
        if 'not' in settings[pclnsig]:
          reportVariant = True

          # Loop over all of the values in the 'not' category. If any of these terms appear in any of the
          # new ClinVar annotations, do not report the variant
          for sig in settings[pclnsig]['not']:
            for clnsig in clnsigs:
              if sig.lower() in clnsig.lower(): reportVariant = False

        # If 'is' occurs in the settings, this means that if the CLNSIG in the new ClinVar vcf 'is'
        # any of these values, the variant should be reported
        if 'is' in settings[pclnsig]:
          for clnsig in clnsigs:
            if clnsig.lower() in [x.lower() for x in settings[pclnsig]['is']]:
              reportVariant = True
              break

    # Store the change for reporting for auditing purposes
    pclnsig = '_'.join(pclnsigs)
    clnsig  = '_'.join(clnsigs)
    if reportVariant:
      if pclnsig not in auditing['true']: auditing['true'][pclnsig] = {}
      if clnsig not in auditing['true'][pclnsig]: auditing['true'][pclnsig][clnsig] = 1
      else: auditing['true'][pclnsig][clnsig] += 1
    else:
      if pclnsig not in auditing['false']: auditing['false'][pclnsig] = {}
      if clnsig not in auditing['false'][pclnsig]: auditing['false'][pclnsig][clnsig] = 1
      else: auditing['false'][pclnsig][clnsig] += 1

    # Convert the CLNREVSTAT fields into the number of stars and then write out in vcf format
    info = 'PCLNSIG=' + str(fields[7]) + ';CLNSIG=' + str(fields[8]) + ';PCLNREVSTAT=' + str(pclnrev) + ';CLNREVSTAT=' + str(clnrev)

    # If the variant is to be reported, write it out
    if reportVariant: print('\t'.join(fields[0:7]), info, sep = '\t', file = outputFile)

  # Close the output file
  outputFile.close()

  # Convert the output vcf file into a compressed, indexed vcf files
  command = bcftools.compressAndIndex(bcftoolsExe, output, str(output) + '.gz', 'z')
  process = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  stdout, stderr = process.communicate()
  process.wait()

  # Write out auditing information
  print()
  print('A total of ', str(totalVariants) + ' variants were processed')
  print()
  print('The following variants were reported:')
  noVariants = 0
  for pclnsig in auditing['true']:
    for clnsig in auditing['true'][pclnsig]:
      print('  ', pclnsig, ' -> ', clnsig, ': ', auditing['true'][pclnsig][clnsig], sep = '') 
      noVariants += auditing['true'][pclnsig][clnsig]
  print('Total variants reported: ', noVariants, sep = '')
  print()
  print('The following variants were NOT reported:')
  noVariants = 0
  for pclnsig in auditing['false']:
    for clnsig in auditing['false'][pclnsig]:
      print('  ', pclnsig, ' -> ', clnsig, ': ', auditing['false'][pclnsig][clnsig], sep = '') 
      noVariants += auditing['false'][pclnsig][clnsig]
  print('Total variants NOT reported: ', noVariants, sep = '')

# Merge the files containing the variants to report, sort, compress and index
def mergeFiles(new, removed, updated, output):
  global workingDir
  global bcftoolsExe

  # Merge the vcf files containing new, removed and updated variants
  inputVcfs = []
  inputVcfs.append(str(workingDir) + str(new))
  inputVcfs.append(str(workingDir) + str(removed))
  inputVcfs.append(str(workingDir) + str(updated))
  command = bcftools.merge(bcftoolsExe, inputVcfs, 'z', str(workingDir) + str(output))
  process = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  stdout, stderr = process.communicate()
  process.wait()

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

# Map the CLNREVSTAT text to a star value
stars = {}
stars['practice_guidline']                                    = "4"
stars['reviewed_by_expert_panel']                             = "3"
stars['criteria_provided,_multiple_submitters,_no_conflicts'] = "2"
stars['criteria_provided,_conflicting_interpretations']       = "1"
stars['criteria_provided,_single_submitter']                  = "1"
stars['no_assertion_criteria_provided']                       = "0"
stars['no_assertion_provided']                                = "0"
stars['no_interpretation_for_the_single_variant']             = "0"

# Keep track of all the CLNSIG updates and whether they were reported
auditing = {}
auditing['true']  = {}
auditing['false'] = {}

if __name__ == "__main__":
  main()
