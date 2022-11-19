#!/usr/bin/python

from __future__ import print_function
from os.path import exists

import os
import math
import argparse
import json
from random import random

from sys import path
path.append("/".join(os.path.dirname(os.path.abspath(__file__)).split("/")[0:-1]) + "/api_commands")
path.append("/".join(os.path.dirname(os.path.abspath(__file__)).split("/")[0:-1]) + "/common_components")
#import mosaic_config
import tools_bcftools as bcf

def main():
  global mosaicConfig

  # Parse the command line
  args = parseCommandLine()

  # Parse the mosaic configuration file
  #mosaicRequired = {"token": True, "url": True, "attributesProjectId": False}
  #mosaicConfig   = mosaic_config.parseConfig(args, mosaicRequired)

  # Check the arguments are valid
  checkArguments(args)

  # Generate intersections of the clinVar files
  generateIntersections(args)

  # Process the outputs of the intersection
 # uniqueToOld()
 # uniqueToNew(args)
  shared(args)

# Input options
def parseCommandLine():
  parser = argparse.ArgumentParser(description='Process the command line arguments')

  # Arguments related to the config file
  #parser.add_argument('--config', '-c', required = False, metavar = "string", help = "A config file containing token / url information")
  #parser.add_argument('--token', '-t', required = False, metavar = "string", help = "The Mosaic authorization token")
  #parser.add_argument('--url', '-u', required = False, metavar = "string", help = "The base url for Mosaic curl commands, up to an including \"api\". Do NOT include a trailing /")

  # The project id to which the filter is to be added is required
  #parser.add_argument('--project_id', '-p', required = True, metavar = "integer", help = "The Mosaic project id to upload attributes to")

  # Arguments required for the ClinVar data
  parser.add_argument('--reference', '-r', required = True, metavar = 'string', help = 'The reference genome to use. Allowed values: ' + ', '.join(allowedReferences))
  parser.add_argument('--data_directory', '-d', required = True, metavar = 'string', help = 'The path to the directory where the resources live')
  parser.add_argument('--old_date', '-o', required = True, metavar = "string", help = 'The date (format: YYYYMMDD) of the old ClinVar file')
  parser.add_argument('--new_date', '-n', required = True, metavar = "string", help = 'The date (format: YYYYMMDD) of the new ClinVar file')

  return parser.parse_args()

# Check the arguments are valid
def checkArguments(args):
  global allowedReferences
  global oldFilename
  global newFilename
  global clinvarPath

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

# Upload the variants from the supplied vcf file
def generateIntersections(args):
  global oldFilename
  global newFilename
  global clinvarPath
  global compDir

  # Generate the name of the output directory. This should be the dates of the old file then the new file separated by _
  compDir = str(args.old_date) + '_' + str(args.new_date)

  # Execute the bcftools command
  try: data = os.popen(bcf.isec(str(clinvarPath) + str(oldFilename), str(clinvarPath) + str(newFilename), compDir)).read()
  except: fail('Unable to execute bcftools isec command')

# Process the vcf file containing variants that are unique to the "old" vcf file. This is variants that have been removed from ClinVar
# between the two supplied time points. In this case, find all pathogenic and VUS variants.
def uniqueToOld():
  global clinvarPath
  global compDir

  # File 0000.vcf.gz contains the variants that are unique to the old vcf file
  for record in os.popen(bcf.query(str(clinvarPath) + str(compDir) + '/0000.vcf.gz', ['CLNSIG'])).readlines():
    fields = record.rstrip().split('\t')

    # Only retain records that have a pathogenic or uncertain variant
    #if 'athogenic' in fields[5] or 'ncertain' in fields[5]:
     # print(fields[5])

# Process the vcf file containing variants that are unique to the "new" vcf file. This is variants that have been added to ClinVar
# between the two supplied time points. In this case, find all pathogenic and VUS variants.
def uniqueToNew(args):
  global clinvarPath
  global compDir
  files = {}

  # File 0001.vcf.gz contains the variants that are unique to the new vcf file
  vcf = str(clinvarPath) + str(compDir) + '/0001.vcf.gz'

  # Get the header from the vcf
  header = os.popen(bcf.getHeader(vcf)).readlines()

  # Open output files to store variants and add the header
  files['pathogenic'] = {'name': str(clinvarPath) + str(compDir) + '/unique_' + str(args.new_date) + '_pathogenic.vcf'}
  files['anyPath']    = {'name': str(clinvarPath) + str(compDir) + '/unique_' + str(args.new_date) + '_other_pathogenic.vcf'}
  files['uncertain']  = {'name': str(clinvarPath) + str(compDir) + '/unique_' + str(args.new_date) + '_uncertain.vcf'}
  files['pathogenic']['file'] = open(files['pathogenic']['name'], 'w')
  files['anyPath']['file']    = open(files['anyPath']['name'], 'w')
  files['uncertain']['file']  = open(files['uncertain']['name'], 'w')
  for line in header:
    for filename in files: print(line.rstrip(), file = files[filename]['file'])

  # Loop over the input vcf file looking for the records of interest
  for record in os.popen(bcf.queryVcf(vcf, ['CLNSIG'])).readlines():
    fields = record.rstrip().split('\t')

    # Only retain records that have a pathogenic or uncertain variant
    if 'Pathogenic' in fields[7]:
      fields[7] = 'CLNSIG=' + str(fields[7])
      print('\t'.join(fields), file = files['pathogenic']['file'])
    elif 'athogenic' in fields[7]:
      fields[7] = 'CLNSIG=' + str(fields[7])
      print('\t'.join(fields), file = files['anyPath']['file'])
    elif 'ncertain' in fields[7]:
      fields[7] = 'CLNSIG=' + str(fields[7])
      print('\t'.join(fields), file = files['uncertain']['file'])

  # Close the output files, the compress and index them and delete the uncompressed files
  finalizeFiles(files)

# Process the vcf file containing variants that are shared by the "old" and "new" vcf files. In this case, we are interested
# in variants that:
#   1. have been upgraded from any non-pathogenic to any pathogenic significance
#   2. have been upgraded from a benign to a VUS
#   3. have been downgraded from any pathognic to VUS or benign significance
def shared(args):
  global clinvarPath
  global compDir

  # Get the name of the "new" ClinVar vcf file
  oldVcf = str(clinvarPath) + str(compDir) + '/0002.vcf.gz'
  newVcf = str(clinvarPath) + str(compDir) + '/0003.vcf.gz'

  # Get the header from the vcf
  header = os.popen(bcf.getHeader(newVcf)).readlines()

  # Open output files to store variants and add the header
  oldFiles = initialiseFiles('old', header)
  newFiles = initialiseFiles('new', header)

  # Process the vcf files to extract files based on the significance
  processVcf(oldVcf, oldFiles)
  processVcf(newVcf, newFiles)

  # Close the output newFiles, the compress and index them and delete the uncompressed newFiles
  finalizeFiles(oldFiles)
  finalizeFiles(newFiles)

  # Now use these vcf files separated on significance to find variants in the desired categories
  os.popen(bcf.isecInt(newFiles['path']['name'] + '.gz', oldFiles['benign']['name'] + '.gz', str(clinvarPath) + str(compDir) + '/benign_path.vcf.gz')).read()
  os.popen(bcf.isecInt(newFiles['path']['name'] + '.gz', oldFiles['conflicting']['name'] + '.gz', str(clinvarPath) + str(compDir) + '/conflicting_path.vcf.gz')).read()
  os.popen(bcf.isecInt(newFiles['path']['name'] + '.gz', oldFiles['vus']['name'] + '.gz', str(clinvarPath) + str(compDir) + '/vus_path.vcf.gz')).read()

  os.popen(bcf.isecInt(newFiles['vus']['name'] + '.gz', oldFiles['benign']['name'] + '.gz', str(clinvarPath) + str(compDir) + '/benign_vus.vcf.gz')).read()
  os.popen(bcf.isecInt(newFiles['conflicting']['name'] + '.gz', oldFiles['benign']['name'] + '.gz', str(clinvarPath) + str(compDir) + '/benign_conflicting.vcf.gz')).read()

  os.popen(bcf.isecInt(newFiles['conflicting']['name'] + '.gz', oldFiles['path']['name'] + '.gz', str(clinvarPath) + str(compDir) + '/path_conflicting.vcf.gz')).read()
  os.popen(bcf.isecInt(newFiles['vus']['name'] + '.gz', oldFiles['path']['name'] + '.gz', str(clinvarPath) + str(compDir) + '/path_vus.vcf.gz')).read()
  os.popen(bcf.isecInt(newFiles['benign']['name'] + '.gz', oldFiles['path']['name'] + '.gz', str(clinvarPath) + str(compDir) + '/path_benign.vcf.gz')).read()

  os.popen(bcf.isecInt(newFiles['vus']['name'] + '.gz', oldFiles['conflicting']['name'] + '.gz', str(clinvarPath) + str(compDir) + '/conflicting_vus.vcf.gz')).read()
  os.popen(bcf.isecInt(newFiles['benign']['name'] + '.gz', oldFiles['conflicting']['name'] + '.gz', str(clinvarPath) + str(compDir) + '/conflicting_benign.vcf.gz')).read()

  os.popen(bcf.isecInt(newFiles['benign']['name'] + '.gz', oldFiles['vus']['name'] + '.gz', str(clinvarPath) + str(compDir) + '/vus_benign.vcf.gz')).read()

# Define and open a set of files for output
def initialiseFiles(name, header):
  files                = {}
  files['conflicting'] = {'name': str(clinvarPath) + str(compDir) + '/' + str(name) + '_conflicting.vcf'}
  files['path']        = {'name': str(clinvarPath) + str(compDir) + '/' + str(name) + '_path.vcf'}
  files['vus']         = {'name': str(clinvarPath) + str(compDir) + '/' + str(name) + '_vus.vcf'}
  files['benign']      = {'name': str(clinvarPath) + str(compDir) + '/' + str(name) + '_benign.vcf'}

  files['conflicting']['file'] = open(files['conflicting']['name'], 'w')
  files['path']['file']        = open(files['path']['name'], 'w')
  files['vus']['file']         = open(files['vus']['name'], 'w')
  files['benign']['file']      = open(files['benign']['name'], 'w')
  for line in header:
    for filename in files: print(line.rstrip(), file = files[filename]['file'])

  # Return the initialised files
  return files

# Loop over the vcf file looking for the records of interest
def processVcf(vcf, files):
  for record in os.popen(bcf.queryVcf(vcf, ['CLNSIG'])).readlines():
    fields = record.rstrip().split('\t')

    # Find variants (from the "new" ClinVar vcf) that are listed as a pathogenic significance. These could potentially
    # have been upgraded.
    if 'Conflicting' in fields[7]:
      fields[7] = 'CLNSIG=' + str(fields[7])
      print('\t'.join(fields), file = files['conflicting']['file'])
    elif 'Pathogenic' in fields[7] or 'Likely_pathogenic' in fields[7]:
      fields[7] = 'CLNSIG=' + str(fields[7])
      print('\t'.join(fields), file = files['path']['file'])
    elif 'Uncertain' in fields[7]:
      fields[7] = 'CLNSIG=' + str(fields[7])
      print('\t'.join(fields), file = files['vus']['file'])
    elif 'Benign' in fields[7] or 'Likely_benign' in fields[7]:
      fields[7] = 'CLNSIG=' + str(fields[7])
      print('\t'.join(fields), file = files['benign']['file'])

# Close output vcf files that have been written to, compress and index them and remove the uncompressed files
def finalizeFiles(files):
  for filename in files:
    files[filename]['file'].close()

    # Compress and index the output files
    try: data = os.popen(bcf.compress(files[filename]['name'], 'z')).read()
    except: fail('Could not compress file ' + str(files[filename]['name']))
    try: data = os.popen(bcf.index(files[filename]['name'] + '.gz')).read()
    except: fail('Could not compress file ' + str(files[filename]['name']))

    # Delete the uncompressed vcf files
    os.remove(files[filename]['name'])

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = '')
  exit(1)

# Initialise global variables
mosaicConfig = {}

# Store info on allowed values
allowedReferences = ['37', '38']

# ClinVar files
clinvarPath = False
compDir     = False
oldFilename = False
newFilename = False

if __name__ == "__main__":
  main()
