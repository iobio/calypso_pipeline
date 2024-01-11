from os.path import exists

import os
import math
import argparse
import json
import subprocess
from random import random

from sys import path
path.append("/".join(os.path.dirname(os.path.abspath(__file__)).split("/")[0:-1]) + "/api_commands")
path.append("/".join(os.path.dirname(os.path.abspath(__file__)).split("/")[0:-1]) + "/common_components")
import tools_bcftools as bcf

def main():

  # Parse the command line
  args = parseCommandLine()

  # Check the arguments are valid
  checkArguments(args)

  # Generate all intersections and perform comparisons to generate the set of vcf files that describe the differences between
  # the previous and the current ClinVar vcf files
  compare(compDir, previousVcf, currentVcf)

# Input options
def parseCommandLine():
  parser = argparse.ArgumentParser(description='Process the command line arguments')

  # Arguments required for the ClinVar data
  parser.add_argument('--reference', '-r', required = True, metavar = 'string', help = 'The reference genome to use. Allowed values: ' + ', '.join(allowedReferences))
  parser.add_argument('--data_directory', '-d', required = True, metavar = 'string', help = 'The path to the directory where the resources live')
  parser.add_argument('--old_date', '-o', required = True, metavar = "string", help = 'The date (format: YYYYMMDD) of the old ClinVar file')
  parser.add_argument('--new_date', '-n', required = True, metavar = "string", help = 'The date (format: YYYYMMDD) of the new ClinVar file')

  return parser.parse_args()

# Check the arguments are valid
def checkArguments(args):
  global allowedReferences

  # Check the reference is allowed
  if args.reference not in allowedReferences:
    message = 'The allowed references (--reference, -r) are:'
    for ref in allowedReferences: message += '\n  ' + str(ref)
    fail(message)

  # Make sure the data directory has a trailing /
  if not args.data_directory.endswith('/'): args.data_directory = args.data_directory + '/'

  # Set the ClinVar vcf file names and check that they exist
  oldFilename = 'clinvar_' + args.old_date + '.vcf.gz'
  newFilename = 'clinvar_' + args.new_date + '.vcf.gz'
  clinvarPath = str(args.data_directory) + 'GRCh' + str(args.reference) + '/clinvar/'
  if not exists(str(clinvarPath) + str(oldFilename)): fail('"Old" vcf file could not be found: ' + str(clinvarPath) + str(oldFilename))
  if not exists(str(clinvarPath) + str(newFilename)): fail('"New" vcf file could not be found: ' + str(clinvarPath) + str(newFilename))

# Generate all intersections and perform comparisons to generate the set of vcf files that describe the differences between
# the previous and the current ClinVar vcf files
def compare(compDir, previousVcf, currentVcf):

  # Open a bash script to write commands to
  bashScript = openBashScript()

  # Generate intersections of the clinVar files
  generateIntersections(compDir, previousVcf, currentVcf, bashScript)

  # Process the outputs of the intersection
  createdFiles = unique(compDir, bashScript)
  createdFiles = shared(compDir, createdFiles, bashScript)

  # Merge all the final files into a single vcf. This will be used to compare with the patient vcf file
  #mergeVcfFiles(createdFiles)

  # Delete all the intermediate files
  deleteIntersectionFiles(compDir, bashScript)

  # Close the bash script and makee executable
  closeBashScript(bashScript)

# Open a bash script to write commands to
def openBashScript():
  bashScript = {'name': 'clinvar_script.sh'}
  bashScript['file'] = open(bashScript['name'], 'w')
  print('#!/bin/bash', file = bashScript['file'])
  print(file = bashScript['file'])

  return bashScript

# Upload the variants from the supplied vcf file
def generateIntersections(compDir, previousVcf, currentVcf, bashScript):

  print('echo -n "Generating intersection and complements between ClinVar files..."', file = bashScript['file'])
  print(bcf.isec(previousVcf, currentVcf, compDir), file = bashScript['file'])
  print('echo "complete"', file = bashScript['file'])

# Process the vcf fileis containing variants that are unique to the vcf files
def unique(compDir, bashScript):
  global header

  # File 0000.vcf.gz contains variants unique to the previous vcf file and 0001.vcf.gz those unique to the current vcf file
  pUniqueVcf = str(compDir) + '/0000.vcf.gz'
  cUniqueVcf = str(compDir) + '/0001.vcf.gz'

  # Get the header from the vcf
  print(file = bashScript['file'])
  print('echo -n "Getting the vcf header..."', file = bashScript['file'])
  header = str(compDir) + '/header.vcf'
  print(bcf.createHeaderVcf(cUniqueVcf, header), file = bashScript['file'])
  print('echo "complete"', file = bashScript['file'])

  # Open output files to store variants and add the header
  pFiles = initialiseFiles(compDir, 'unique_previous', header, bashScript)
  cFiles = initialiseFiles(compDir, 'unique_current', header, bashScript)

  # Process the vcf files to extract files based on the significance
  print(file = bashScript['file'])
  print('echo -n "Processing unique complements of the ClinVar vcf files..."', file = bashScript['file'])
  processVcf(pUniqueVcf, pFiles, bashScript)
  processVcf(cUniqueVcf, cFiles, bashScript)
  print('echo "complete"', file = bashScript['file'])

  # Compress and index the output files and delete the uncompressed files
  print(file = bashScript['file'])
  print('echo -n "Compressing and indexing unique complements..."', file = bashScript['file'])
  finalizeFiles(pFiles, bashScript)
  finalizeFiles(cFiles, bashScript)
  print('echo "complete"', file = bashScript['file'])

  # Keep track of all created files
  createdFiles = []
  for pFile in pFiles: createdFiles.append(pFiles[pFile]['name'])
  for cFile in cFiles: createdFiles.append(cFiles[cFile]['name'])

  # Return the list of created files
  return createdFiles

# Process the vcf file containing variants that are shared by the "previous (p prefix)" and "current (c prefix)" vcf files.
# In this case, we are interested in variants that:
#
#   1. have been upgraded from any non-pathogenic to any pathogenic significance
#   2. have been upgraded from a benign to a VUS
#   3. have been downgraded from any pathognic to VUS or benign significance
def shared(compDir, createdFiles, bashScript):
  global header

  # Get the name of the "new" ClinVar vcf file
  pSharedVcf = str(compDir) + '/0002.vcf.gz'
  cSharedVcf = str(compDir) + '/0003.vcf.gz'

  # Open output files to store variants and add the header
  pFiles = initialiseFiles(compDir, 'previous', header, bashScript)
  cFiles = initialiseFiles(compDir, 'current', header, bashScript)

  # Process the vcf files to extract files based on the significance
  print(file = bashScript['file'])
  print('echo -n "Processing intersection of the ClinVar vcf files..."', file = bashScript['file'])
  processVcf(pSharedVcf, pFiles, bashScript)
  processVcf(cSharedVcf, cFiles, bashScript)
  print('echo "complete"', file = bashScript['file'])

  # Close the files, the compress and index them and delete the uncompressed files
  print(file = bashScript['file'])
  print('echo -n "Compressing and indexing intersections..."', file = bashScript['file'])
  finalizeFiles(pFiles, bashScript)
  finalizeFiles(cFiles, bashScript)
  print('echo "complete"', file = bashScript['file'])

  # Now use these vcf files separated on significance to find variants in the desired categories and index the created files
  print(file = bashScript['file'])
  print('echo -n "Determining ClinVar updates..."', file = bashScript['file'])
  createdFiles = getUpdateVcf(pFiles['benign']['name'] + '.gz', cFiles['path']['name'] + '.gz', str(compDir) + '/benign_path.vcf.gz', createdFiles, bashScript)
  createdFiles = getUpdateVcf(pFiles['conflicting']['name'] + '.gz', cFiles['path']['name'] + '.gz', str(compDir) + '/conflicting_path.vcf.gz', createdFiles, bashScript)
  createdFiles = getUpdateVcf(pFiles['vus']['name'] + '.gz', cFiles['path']['name'] + '.gz', str(compDir) + '/vus_path.vcf.gz', createdFiles, bashScript)

  createdFiles = getUpdateVcf(pFiles['benign']['name'] + '.gz', cFiles['vus']['name'] + '.gz', str(compDir) + '/benign_vus.vcf.gz', createdFiles, bashScript)
  createdFiles = getUpdateVcf(pFiles['conflicting']['name'] + '.gz', cFiles['vus']['name'] + '.gz', str(compDir) + '/conflicting_vus.vcf.gz', createdFiles, bashScript)
  createdFiles = getUpdateVcf(pFiles['path']['name'] + '.gz', cFiles['vus']['name'] + '.gz', str(compDir) + '/path_vus.vcf.gz', createdFiles, bashScript)

  createdFiles = getUpdateVcf(pFiles['benign']['name'] + '.gz', cFiles['conflicting']['name'] + '.gz', str(compDir) + '/benign_conflicting.vcf.gz', createdFiles, bashScript)
  createdFiles = getUpdateVcf(pFiles['vus']['name'] + '.gz', cFiles['conflicting']['name'] + '.gz', str(compDir) + '/vus_conflicting.vcf.gz', createdFiles, bashScript)
  createdFiles = getUpdateVcf(pFiles['path']['name'] + '.gz', cFiles['conflicting']['name'] + '.gz', str(compDir) + '/path_conflicting.vcf.gz', createdFiles, bashScript)

  createdFiles = getUpdateVcf(pFiles['benign']['name'] + '.gz', cFiles['conflicting']['name'] + '.gz', str(compDir) + '/benign_conflicting.vcf.gz', createdFiles, bashScript)
  createdFiles = getUpdateVcf(pFiles['vus']['name'] + '.gz', cFiles['conflicting']['name'] + '.gz', str(compDir) + '/vus_conflicting.vcf.gz', createdFiles, bashScript)
  createdFiles = getUpdateVcf(pFiles['path']['name'] + '.gz', cFiles['conflicting']['name'] + '.gz', str(compDir) + '/path_conflicting.vcf.gz', createdFiles, bashScript)
  print('echo "complete"', file = bashScript['file'])

  # Remove the intermediate files created in this step
  removeFiles(pFiles, bashScript)
  removeFiles(cFiles, bashScript)

  # Return the updated list of create files
  return createdFiles

# Define and open a set of files for output
def initialiseFiles(compDir, name, header, bashScript):
  files                = {}
  files['conflicting'] = {'name': str(compDir) + '/' + str(name) + '_conflicting.vcf'}
  files['path']        = {'name': str(compDir) + '/' + str(name) + '_path.vcf'}
  files['vus']         = {'name': str(compDir) + '/' + str(name) + '_vus.vcf'}
  files['benign']      = {'name': str(compDir) + '/' + str(name) + '_benign.vcf'}

  #files['conflicting']['file'] = open(files['conflicting']['name'], 'w')
  #files['path']['file']        = open(files['path']['name'], 'w')
  #files['vus']['file']         = open(files['vus']['name'], 'w')
  #files['benign']['file']      = open(files['benign']['name'], 'w')

  for filename in files: print('cp ', header, ' ', files[filename]['name'], sep = '', file = bashScript['file'])

  # Return the initialised files
  return files

# Loop over the vcf file looking for the records of interest
def processVcf(vcf, files, bashScript):

  # Get the bcftool query command to include only the CLNSIGN info
  print(bcf.queryVcf(vcf, ['CLNSIG']), ' \\', sep = '', file = bashScript['file'])
  print('  | awk \'BEGIN {FS="\\t"; OFS="\\t"} { \\', file = bashScript['file'])
  print('  if ($8 ~ /Conflicting/) {print $1,$2,$3,$4,$5,$6,$7,"CLNSIG="$8 >> "', files['conflicting']['name'], '";} \\', sep = '', file = bashScript['file'])
  print('  else if ($8 ~ /athogenic/) {print $1,$2,$3,$4,$5,$6,$7,"CLNSIG="$8 >> "', files['path']['name'], '";} \\', sep = '', file = bashScript['file'])
  print('  else if ($8 ~ /Uncertain/) {print $1,$2,$3,$4,$5,$6,$7,"CLNSIG="$8 >> "', files['vus']['name'], '";} \\', sep = '', file = bashScript['file'])
  print('  else if ($8 ~ /enign/) {print $1,$2,$3,$4,$5,$6,$7,"CLNSIG="$8 >> "', files['benign']['name'], '";} \\', sep = '', file = bashScript['file'])
  print('  }\'', file = bashScript['file'])

# Close output vcf files that have been written to, compress and index them and remove the uncompressed files
def finalizeFiles(files, bashScript):
  for filename in files:

    # Compress and index the output files
    print(bcf.compress(files[filename]['name'], 'z'), file = bashScript['file'])
    print(bcf.index(files[filename]['name'] + '.gz'), file = bashScript['file'])
    print('rm -f ', files[filename]['name'], sep = '', file = bashScript['file'])
    print(file = bashScript['file'])

# Generate a vcf containing all variants that were updated from one significance to another
def getUpdateVcf(pFile, cFile, uFile, createdFiles, bashScript):

  # Use subprocess.call to wait until the intersection is complete before indexing
  print(bcf.isecInt(cFile, pFile, uFile), file = bashScript['file'])
  print(bcf.index(uFile), file = bashScript['file'])

  # Return an updated list of created files
  createdFiles.append(uFile)
  return createdFiles

# Remove created files
def removeFiles(files, bashScript):
  for filename in files:
    print('rm -f ', files[filename]['name'], '.gz', sep = '', file = bashScript['file'])
    print('rm -f ', files[filename]['name'], '.gz.tbi', sep = '', file = bashScript['file'])

# Merge all the final files into a single vcf. This will be used to compare with the patient vcf file
def mergeVcfFiles(createdFiles):
  print(bcf.merge(createdFiles, 'z', str(compDir) + '/merged.vcf.gz'))
  os.popen(bcf.merge(createdFiles, 'z', str(compDir) + '/merged.vcf.gz'))
  os.popen(bcf.index(str(compDir) + '/merged.vcf.gz'))

# Now delete all the intermediate files
def deleteIntersectionFiles(compDir, bashScript):
  global header

  print('rm -f ', str(compDir), '/0000.vcf.gz', sep = '', file = bashScript['file'])
  print('rm -f ', str(compDir), '/0000.vcf.gz.tbi', sep = '', file = bashScript['file'])
  print('rm -f ', str(compDir), '/0001.vcf.gz', sep = '', file = bashScript['file'])
  print('rm -f ', str(compDir), '/0001.vcf.gz.tbi', sep = '', file = bashScript['file'])
  print('rm -f ', str(compDir), '/0002.vcf.gz', sep = '', file = bashScript['file'])
  print('rm -f ', str(compDir), '/0002.vcf.gz.tbi', sep = '', file = bashScript['file'])
  print('rm -f ', str(compDir), '/0003.vcf.gz', sep = '', file = bashScript['file'])
  print('rm -f ', str(compDir), '/0003.vcf.gz.tbi', sep = '', file = bashScript['file'])
  print('rm -f ', str(compDir), '/README.txt', sep = '', file = bashScript['file'])
  print('rm -f ', str(compDir), '/sites.txt', sep = '', file = bashScript['file'])
  print('rm -f ', header, sep = '', file = bashScript['file'])
  print('rm -f ', bashScript['name'], sep = '', file = bashScript['file'])

# Close the bash script
def closeBashScript(bashScript):
  bashScript['file'].close()
  makeExecutable = os.popen('chmod +x ' + bashScript['name']).read()

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = '')
  exit(1)

# Store info on allowed values
allowedReferences = ['37', '38']

# Store the name of the header file
header = False

if __name__ == "__main__":
  main()
