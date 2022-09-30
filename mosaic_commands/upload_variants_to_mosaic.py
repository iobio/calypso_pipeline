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
import api_variants as api_v

def main():
  global mosaicConfig

  # Parse the command line
  args = parseCommandLine()

  # Parse the Mosaic config file to get the token and url for the api calls
  mosaicRequired = {"token": True, "url": True, "attributesProjectId": False}
  mosaicConfig   = mosaic_config.parseConfig(args, mosaicRequired)

  # Upload variants
  uploadVariants(args)

# Input options
def parseCommandLine():
  global version

  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  parser.add_argument('--vcf_file', '-i', required = True, metavar = "string", help = "The vcf file containing variants to upload")
  parser.add_argument('--project', '-p', required = True, metavar = "integer", help = "The Mosaic project id to upload attributes to")
  parser.add_argument('--description', '-d', required = True, metavar = "string", help = "A description of the variants being uploaded")

  # Optional mosaic arguments
  parser.add_argument('--config', '-c', required = False, metavar = "string", help = "The config file for Mosaic")
  parser.add_argument('--token', '-t', required = False, metavar = "string", help = "The Mosaic authorization token")
  parser.add_argument('--url', '-u', required = False, metavar = "string", help = "The base url for Mosaic")

  # Version
  parser.add_argument('--version', '-v', action="version", version="Calypso variant uploader version: " + version)

  return parser.parse_args()

# Upload the supplied VCF file to Mosaic
def uploadVariants(args):
  global mosaicConfig

  try: data = os.popen(api_v.postUploadVariants(mosaicConfig, args.vcf_file, "allele", args.description, args.project))
  except: fail("Unable to upload file")

# Initialise global variables. These annotations are in the order they should be output to file
version = "0.0.2"
mosaicConfig = {}

if __name__ == "__main__":
  main()
