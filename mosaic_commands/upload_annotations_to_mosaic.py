#!/usr/bin/python

from __future__ import print_function
from os.path import exists

import os
import argparse
import json

def main():

  # Parse the command line
  args = parseCommandLine()

  # Parse the config file
  parseConfig(args)

  # Upload variants
  uploadAnnotations(args)

# Input options
def parseCommandLine():
  global version

  parser = argparse.ArgumentParser(description='Process the command line')
  parser.add_argument('--config', '-c', required = True, metavar = "string", help = "A config file containing token / url information")
  parser.add_argument('--api_commands', '-a', required = True, metavar = "string", help = "The path to the directory of api commands")
  parser.add_argument('--tsv_file', '-i', required = True, metavar = "string", help = "The vcf file containing variants to upload")
  parser.add_argument('--project', '-p', required = True, metavar = "integer", help = "The Mosaic project id to upload attributes to")
  parser.add_argument('--version', '-v', action="version", version="Calypso variant uploader version: " + version)

  return parser.parse_args()

# Parse the config file
def parseConfig(args):
  global token
  global apiUrl

  # Check that the config file exists
  if not exists(args.config):
    print("Config file does not exist")
    exit(1)

  # Parse the file and extract recognised fields
  with open(args.config) as configFile:
    for line in configFile:
      if "=" in line:
        argument = line.rstrip().split("=")

        # Strip any whitespace from the arguments
        argument[0] = argument[0].strip()
        argument[1] = argument[1].strip()

        # Set the recognized values
        if argument[0] == "MOSAIC_TOKEN": token = argument[1]
        if argument[0] == "MOSAIC_URL":
          if argument[1].endswith("/"): apiUrl = argument[1] + "api"
          else: apiUrl = argument[1] + "/api"

# Upload the supplied VCF file to Mosaic
def uploadAnnotations(args):
  global token
  global apiUrl

  outFile  = "calypso_upload_annotations.stdout"
  errFile  = "calypso_upload_annotations.stderr"
  command  = args.api_commands + "/upload_variant_annotations.sh " + str(token) + " \"" + str(apiUrl) + "\" \"" + str(args.project) + "\" \""
  command += str(args.tsv_file) + "\" > " + str(outFile) + " 2> " + str(errFile)
  data     = os.popen(command).read()

# Initialise global variables. These annotations are in the order they should be output to file
version = "0.0.1"
token   = False
apiUrl  = False

if __name__ == "__main__":
  main()
