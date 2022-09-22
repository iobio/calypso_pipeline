#!/usr/bin/python

from __future__ import print_function
from os.path import exists

import os
import argparse
import json

# Add the path of the common functions and import them
from sys import path
path.append("/".join(os.path.dirname(os.path.abspath(__file__)).split("/")[0:-1]) + "/mosaic_commands")
import mosaic_config
import api_variant_annotations as api_va

def main():
  global mosaicConfig

  # Parse the command line
  args = parseCommandLine()

  # Parse the Mosaic config file to get the token and url for the api calls
  mosaicRequired = {"token": True, "url": True, "attributesProjectId": False}
  mosaicConfig   = mosaic_config.parseConfig(args, mosaicRequired)

  # Upload variants
  uploadAnnotations(args)

# Input options
def parseCommandLine():
  global version

  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  parser.add_argument('--tsv_file', '-i', required = True, metavar = "string", help = "The vcf file containing variants to upload")
  parser.add_argument('--project', '-p', required = True, metavar = "integer", help = "The Mosaic project id to upload attributes to")

  # Optional arguments
  parser.add_argument('--no_deletion', '-d', required = False, action = "store_true", help = "If set, blank values will NOT overwite existing annotation values")

  # Optional mosaic arguments
  parser.add_argument('--config', '-c', required = False, metavar = "string", help = "The config file for Mosaic")
  parser.add_argument('--token', '-t', required = False, metavar = "string", help = "The Mosaic authorization token")
  parser.add_argument('--url', '-u', required = False, metavar = "string", help = "The base url for Mosaic")

  # Version
  parser.add_argument('--version', '-v', action="version", version="Calypso variant uploader version: " + version)

  return parser.parse_args()

# Upload the supplied VCF file to Mosaic
def uploadAnnotations(args):
  global mosaicConfig

  # By default overwrite existing annotations with a blank
  allowDeletion = "false" if args.no_deletion else "true"

  try: data = os.popen(api_va.postUploadVariantAnnotations(mosaicConfig, args.tsv_file, allowDeletion, args.project)).read()
  except: fail("Unable to upload file")

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)

# Initialise global variables. These annotations are in the order they should be output to file
version = "0.0.2"
token   = False
apiUrl  = False
mosaicConfig = {}

if __name__ == "__main__":
  main()
