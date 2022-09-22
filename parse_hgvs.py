#!/usr/bin/python

from __future__ import print_function
from datetime import date
from os.path import exists

import os
import argparse
import json
import math

def main():

  # Parse the command line
  args = parseCommandLine()

  # Parse the Mosaic json file describing the mosaic information for uploading annotations
  parseMosaicJson(args)

  # Parse the HGVS file and update
  parseHgvs(args)

  # Write out any warnings
  writeWarnings(args)

# Input options
def parseCommandLine():
  global version

  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  parser.add_argument('--input_tsv', '-i', required = True, metavar = "string", help = "The input HGVS tsv file")
  parser.add_argument('--reference', '-r', required = True, metavar = "string", help = "The reference genome to use. Allowed values: '37', '38'")

  # Optional mosaic arguments
  parser.add_argument('--mosaic_json', '-m', required = False, metavar = "string", help = "The json file describing the Mosaic parameters")
  parser.add_argument('--data_directory', '-d', required = False, metavar = "string", help = "The path to the directory where the resources live")

  # Version
  parser.add_argument('--version', '-v', action="version", version='SpliceAI parser version: ' + str(version))

  return parser.parse_args()

# Parse the Mosaic json file describing the mosaic information for uploading annotations
def parseMosaicJson(args):
  global mosaicInfo

  # Either the data directory or the mosaic json need to be defined.
  if not args.data_directory and not args.mosaic_json: fail("Either the Mosaic json file (--mosaic_json, -m), or the data directory (--data_directory, -d) must be specified")

  # Define the name of the Mosaic json file. Use the file provided on the command
  # line, and resort to the default file if not included
  if args.mosaic_json: mosaicFilename = args.mosaic_json
  else:

    # Ensure the data path ends with a "/", then add the refrence directory
    if args.data_directory[-1] != "/": args.data_directory += "/"

    if args.reference == '37': mosaicFilename = args.data_directory + "resources_mosaic_GRCh37.json"
    elif args.reference == '38': mosaicFilename = args.data_directory + "resources_mosaic_GRCh38.json"

  # Try and open the file
  try: mosaicFile = open(mosaicFilename, "r")
  except: fail("The file describing Mosaic related information (" + str(mosaicFilename) + ") could not be found")

  # Extract the json information
  try: mosaicData = json.loads(mosaicFile.read())
  except: fail("The json file (" + str(mosaicFilename) + ") is not valid")

  # Store the data version
  try: mosaicInfo["version"] = mosaicData['version']
  except: fail("The Mosaic json (" + str(mosaicFilename) + ") does not include a version")

  # Check that the resource json reference matches the selected reference
  try: mosaicReference = mosaicData['reference']
  except: fail("The Mosaic json does not include a reference genome")
  isRefMatch = False
  if args.reference == '37' and mosaicReference == '37': isRefMatch = True
  elif args.reference == '38' and mosaicReference == '38': isRefMatch = True
  if not isRefMatch: fail("The selected reference (" + str(args.reference) + ") does not match the Mosaic json reference (" + str(mosaicReference) + ")")

  # Get the resources
  try: resources = mosaicData["resources"]
  except: fail("The Mosaic json (" + str(mosaicFilename) + ") does not include a list of resources")
  mosaicInfo["resources"] = {}
  for resource in resources:

    # If a resource is duplicated, throw an error
    if resource in mosaicInfo["resources"]: fail(str(resource) + " appears multiple times in the Mosaic json. Ensure each resource appears only once.")
    mosaicInfo["resources"][resource] = {}

    try: annotation_type = resources[resource]["annotation_type"]
    except: fail("The Mosaic json file does not include the 'annotation_type' for resource: " + str(resource))
    mosaicInfo["resources"][resource]["annotation_type"] = annotation_type
    if annotation_type != "private" and annotation_type != "public": fail("Annotation_type for " + str(resource) + " must be 'public' or 'private'")

    # Check for the "project_id". This is required for public annotations
    try: project_id = resources[resource]["project_id"]
    except: project_id = False
    mosaicInfo["resources"][resource]["project_id"] = project_id
    if annotation_type == "public" and not project_id: fail("The Mosaic json file does not include a 'project_id' for public resource: " + str(resource))

    # Check if there is a delimiter present. This is present if the annotation contains multiple values and needs to
    # be split to extract each annotation (e.g. SpliceAI=ALLELE|SYMBOL|DS_AL|..., would require the delimiter to
    # be set to "|" and each annotation must include a "position". The "INFO" field must also be present to describe
    # the INFO field value that will be split, since the annotation name describes each annotation, but this will
    # not appear in the INFO.
    hasDelimiter = False
    hasInfo      = False
    if "delimiter" in resources[resource]:
      mosaicInfo["resources"][resource]["delimiter"] = resources[resource]["delimiter"]
      hasDelimiter = True
    else: mosaicInfo["resources"][resource]["delimiter"] = False
    if "INFO" in resources[resource]:
      mosaicInfo["resources"][resource]["info"] = resources[resource]["INFO"]
      hasInfo = True
    else: mosaicInfo["resources"][resource]["info"] = False

    if hasDelimiter and not hasInfo: fail("'" + str(resource) + "' has 'delimiter' set, but no 'INFO'. Both or neither of these need to be set")
    if hasInfo and not hasDelimiter: fail("'" + str(resource) + "' has 'INFO' set, but no 'delimiter'. Both or neither of these need to be set")

    # Loop over the annotation information
    mosaicInfo["resources"][resource]["annotations"] = {}
    for annotation in resources[resource]["annotations"]:

      # Read in the required fields
      try: uid = resources[resource]["annotations"][annotation]["uid"]
      except: fail("The Mosaic json does not contain the 'uid' field for annotation '" + str(annotation) + "' for resource '" + str(resource) + "'")
      try: annType = resources[resource]["annotations"][annotation]["type"]
      except: fail("The Mosaic json does not contain the 'type' field for annotation '" + str(annotation) + "' for resource '" + str(resource) + "'")

      # Some annotations require an operation to be performed, e.g. find the maximum of the annotations. If this is the case
      # the check on required fields due to the presence of a delimiter can be ignored, since this annotation is not pulling
      # from the INFO field, but is constructing it
      if "operation" in resources[resource]["annotations"][annotation]:
        operation = resources[resource]["annotations"][annotation]["operation"]

        # If the operation to be performed is finding the max of certain values, a list of positions to find the max of must be supplied
        if operation == "max":
          try: maxPositions = resources[resource]["annotations"][annotation]["positions"]
          except: fail("Annotation '" + str(annotation) + "' for '" + str(resource) + "' requires the \"positions\" list defining the fields to find the max of")

          # Check the positions is a list
          if isinstance(maxPositions, list):
            mosaicInfo["resources"][resource]["annotations"][annotation] = {"uid": uid, "type": annType, "operation": operation, "positions": maxPositions}
          else: fail("Annotation '" + str(annotation) + "' for '" + str(resource) + "' include \"positions\" defining the fields to find the max of, but this needs to be a list")

      # If the delimiter is set, a "position" field must be present with a numerical value.
      else:
        if mosaicInfo["resources"][resource]["delimiter"]:
          try: position = resources[resource]["annotations"][annotation]["position"]
          except: fail("Annotation '" + str(annotation) + "' for '" + str(resource) + "' requires the \"position\" field as the delimiter is set")
          if not isinstance(position, int): fail("Annotation '" + annotation + "' for '" + str(resource) + "' has a position field that is not an integer")
          mosaicInfo["resources"][resource]["annotations"][annotation] = {"uid": uid, "type": annType, "position": position}
        else: mosaicInfo["resources"][resource]["annotations"][annotation] = {"uid": uid, "type": annType}

# Parse the HGVS tsv and update ready for upload to Mosaic
def parseHgvs(args):
  global mosaicInfo

  # Open the tsv file
  try: hgvsTsv = open(args.input_tsv, "r")
  except: fail("The file (" + str(args.input_tsv) + ") does not exist")

  # Open the output file
  outFilename = args.input_tsv + ".tmp"
  if "/" in outFilename: outFilename = outFilename.split("/")[-1]
  try: outFile = open(outFilename, "w")
  except: fail("The file (" + str(outFilename) + ") cannot be opened")

  # Loop over each entry in the file
  for line in hgvsTsv:
    line = line.rstrip()

    # Write the header line back out
    if line.startswith("CHROM"): print(line, file = outFile)
    else:
      fields = line.split("\t")

      # Store the csq field
      csq = fields[5]

      # Check that the chromosome does not include "chr" prefix
      outputLine = fields[0][3:] if fields[0].startswith("chr") else fields[0]

      # Add one to the end position
      outputLine += "\t" + str(fields[1]) + "\t" + str(int(fields[2]) + 1) + "\t" + fields[3] + "\t" + fields[4]

      # Check if the line has any values. If every field is a '.', then skip this line
      hasValue = False
      for index in range(5, len(fields)):
        if fields[index] != ".":
          hasValue = True
          break

      # Only continue if at least one of the annotations has a value
      if hasValue:

        # Get all the available CSQ values
        if ',' in csq: csqValues = csq.split(",")
        else: csqValues = [csq]
        annotationLine = ""

        # Loop over the HGVS annotations
        for annotation in mosaicInfo["resources"]["HGVS"]["annotations"]:
          annotationInfo   = mosaicInfo["resources"]["HGVS"]["annotations"][annotation]
          outputAnnotation = {}

          # Get the position of the annotation in the CSQ field
          position  = annotationInfo["position"]

          # Loop over the availableValues and build up a string of all values for this annotation
          values = []
          for csqValue in csqValues:
            value = csqValue.split("|")[int(position) - 1]

#########
#########
######### REMOVE values that are longer than 255 characters. Need a better solution
#########
#########
            if len(value) > 255: value = ""
            if value: values.append(value)
          annotationLine += "\t" + str(",".join(values))

        # Write out the line
        print(outputLine, annotationLine, sep = "", file = outFile)

  # Close the files
  hgvsTsv.close()
  outFile.close()

  # Overwrite the original tsv with the new file
  command = "mv " + outFilename + " " + args.input_tsv
  data    = os.popen(command).read()

# Write out any warnings
def writeWarnings(args):
  global warnings

  if len(warnings) > 0:
    for warning in warnings:
      print("### WARNINGS ###")
      print()
      print("## ", warning, sep = "")
      print("  ", warnings[warning]["description"], sep = "")
      for value in warnings[warning]["messages"]: print("    ", value, sep = "")

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)

# Initialise global variables
mosaicInfo = {}

# Pipeline version
version = "0.1.1"

# Store warnings to be output at the end
warnings = {}

if __name__ == "__main__":
  main()
