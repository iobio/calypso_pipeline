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

  # Check the supplied parameters are as expected
  checkArguments(args)

  # Check all resources
  checkResources(args)

  # Parse the Mosaic json
  parseMosaicJson(args)

  # Build the toml file for vcfanno
  buildToml()

  # Generate the bash script to run the annotation pipeline
  finalVcf = genBashScript(args)

  # Parse the Mosaic config file to get the token and url for the api calls
  parseConfig(args)

  # Generate script to upload filtered variants to Mosaic
  uploadVariants(args)

  # Generate scripts to upload annotations
  uploadAnnotations(args)

  # Output summary file
  calypsoSummary(args, finalVcf)

  # Write out any warnings
  writeWarnings(args)

# Input options
def parseCommandLine():
  global version

  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  parser.add_argument('--reference', '-r', required = True, metavar = "string", help = "The reference genome to use. Allowed values: '37', '38'")
  parser.add_argument('--family_type', '-f', required = True, metavar = "string", help = "The familty structure. Allowed values: 'singleton', 'duo', 'trio'")
  parser.add_argument('--data_directory', '-d', required = True, metavar = "string", help = "The path to the directory where the resources live")
  parser.add_argument('--input_vcf', '-i', required = True, metavar = "string", help = "The input vcf file to annotate")

  # Optional pipeline arguments
  parser.add_argument('--ped', '-p', required = False, metavar = "string", help = "The pedigree file for the family. Not required for singletons")
  parser.add_argument('--resource_json', '-j', required = False, metavar = "string", help = "The json file describing the annotation resources")

  # Optional mosaic arguments
  parser.add_argument('--mosaic_json', '-m', required = False, metavar = "string", help = "The json file describing the Mosaic parameters")
  parser.add_argument('--config', '-c', required = False, metavar = "string", help = "The config file for Mosaic")
  parser.add_argument('--project_id', '-t', required = False, metavar = "string", help = "The project id that variants will be uploaded to")

  # Version
  parser.add_argument('--version', '-v', action="version", version='Calypso annotation pipline version: ' + str(version))

  return parser.parse_args()

# Check the supplied arguments
def checkArguments(args):
  global allowedReferences
  global allowedFamily
  global isFamily

  # Check the reference is allowed
  if args.reference not in allowedReferences:
    print("The allowed references (--reference, -r) are:")
    for ref in allowedReferences: print("  ", ref, sep = "")
    exit(1)

  # Check the family type is allowed
  if args.family_type not in allowedFamily:
    print("The allowed family types (--family-type, -f) are:")
    for family in allowedFamily: print("  ", family, sep = "")
    exit(1)

  # If this is a singleton, no family based analyis should be performed
  if args.family_type == 'singleton': isFamily = False

  # Check that the input files exist and have the correct extensions
  if not exists(args.input_vcf): fail("The vcf file could not be found")
  elif not args.input_vcf.endswith(".vcf.gz"): fail("The input vcf file (--vcf, -v) must be a compressed, indexed vcf and have the extension '.vcf.gz'")
  if isFamily:
    if not args.ped: fail("A ped file needs to specified (--ped, -p) for family type \"" + str(args.family_type) + "\"")
    if not exists(args.ped): fail("The ped file could not be found")
    elif not args.ped.endswith(".ped"): fail("The input ped file (--ped, -p) must have the extension '.ped'")

# Parse the json file describing the resources for the selected genome build,
# check the files exist, and store the versions.
def checkResources(args):
  global resourceInfo

  # Ensure the data path ends with a "/", then add the refrence directory
  if args.data_directory[-1] != "/": args.data_directory += "/"
  resourceInfo["path"] = args.data_directory + "GRCh" + str(args.reference) + "/"

  # Define the name of the resource description json file. Use the file provided on the command
  # line, and resort to the default file if not included
  if args.resource_json: resourceFilename = args.resource_json
  else:
    if args.reference == '37': resourceFilename = args.data_directory + "resources_GRCh37.json"
    elif args.reference == '38': resourceFilename = args.data_directory + "resources_GRCh38.json"

  # Try and open the file
  try: resourceFile = open(resourceFilename, "r")
  except: fail("The file describing the resource files (" + str(resourceFilename) + ") could not be found")

  # Extract the json information
  try: resourceData = json.loads(resourceFile.read())
  except: fail("The json file (" + str(resourceFilename) + ") is not valid")

  # Store the data version
  try: resourceInfo["version"] = resourceData['version']
  except: fail("The resource json file (" + str(resourceFilename) + ") does not include a version")

  # Check that the resource json reference matches the selected reference
  try: resourceReference = resourceData['reference']
  except: fail("The resource json does not include a reference genome")
  isRefMatch = False
  if args.reference == '37' and resourceReference == '37': isRefMatch = True
  elif args.reference == '38' and resourceReference == '38': isRefMatch = True
  if not isRefMatch: fail("The selected reference (" + str(args.reference) + ") does not match the resource json reference (" + str(resourceReference) + ")")

  # Get the resources
  try: resources = resourceData["resources"]
  except: fail("The resources json does not include a list of resources")
  resourceInfo["resources"] = {}
  for resource in resources:

    # If a resource is duplicated, throw an error
    if resource in resourceInfo["resources"]: fail(str(resource) + " appears multiple times in the resources json. Ensure each resource appears only once.")
    resourceInfo["resources"][resource] = {}

    # Get required information for each resource
    try: resourceInfo["resources"][resource]["version"] = resources[resource]["version"]
    except: fail("Version for resource \"" + str(resource) + "\" was not included in the resources json")
    try: resourceInfo["resources"][resource]["file"] = resourceInfo["path"] + resources[resource]["file"]
    except: fail("File for resource \"" + str(resource) + "\" was not included in the resources json")
    try: resourceInfo["resources"][resource]["toml"] = resources[resource]["toml"]
    except: fail("File for resource \"" + str(resource) + "\" was not included in the resources json")
    try: resourceInfo["resources"][resource]["upload"] = resources[resource]["upload_to_mosaic"]
    except: fail("The resources json did not indicate if resource \"" + str(resource) + "\" should be uploaded to Mosaic")

    # If the resource is to be included in a toml file, the fields in the vcf INFO need to be specified
    if resourceInfo["resources"][resource]["toml"]:
      if "fields" in resources[resource]: resourceInfo["resources"][resource]["fields"] = resources[resource]["fields"]
      if "columns" in resources[resource]: resourceInfo["resources"][resource]["columns"] = resources[resource]["columns"]
      if "names" in resources[resource]: resourceInfo["resources"][resource]["names"] = resources[resource]["names"]

    # Check that the file exists
    if not exists(resourceInfo["resources"][resource]["file"]): fail("Resource file " + str(resourceInfo["resources"][resource]["file"]) + " does not exist")

# Parse the Mosaic json file describing the mosaic information for uploading annotations
def parseMosaicJson(args):
  global mosaicInfo

  # Ensure the data path ends with a "/", then add the refrence directory
  if args.data_directory[-1] != "/": args.data_directory += "/"

  # Define the name of the Mosaic json file. Use the file provided on the command
  # line, and resort to the default file if not included
  if args.mosaic_json: mosaicFilename = args.mosaic_json
  else:
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

    # Loop over the annotation information
    mosaicInfo["resources"][resource]["annotations"] = {}
    for annotation in resources[resource]["annotations"]:

      # Read in the required fields
      try: uid = resources[resource]["annotations"][annotation]["uid"]
      except: fail("The Mosaic json does not contain the 'uid' field for annotation '" + str(annotation) + "' for resource '" + str(resource) + "'")
      try: annType = resources[resource]["annotations"][annotation]["type"]
      except: fail("The Mosaic json does not contain the 'type' field for annotation '" + str(annotation) + "' for resource '" + str(resource) + "'")
      mosaicInfo["resources"][resource]["annotations"][annotation] = {"uid": uid, "type": annType}

# Build the toml file
def buildToml():
  global resourceInfo

  # Create a toml file
  try: tomlFile = open("calypso_annotations.toml", "w")
  except: fail("There was a problem opening a file (calypso_annotations.toml) to write to")

  # Add each required resource to the toml file
  for resource in resourceInfo["resources"]:
    if resourceInfo["resources"][resource]["toml"]:
      print("[[annotation]]", file = tomlFile)
      print("file=\"", resourceInfo["resources"][resource]["file"], "\"", sep = "", file = tomlFile)

      # Some of the annotations include "fields", "columns", and "names". Include these in the toml if they exist
      if "fields" in resourceInfo["resources"][resource]:
        text, noValues = tomlInfo(resource, "fields")
        print(text, file = tomlFile)
      if "columns" in resourceInfo["resources"][resource]:
        text, noValues = tomlInfo(resource, "columns")
        print(text, file = tomlFile)
      if "names" in resourceInfo["resources"][resource]:
        text, noValues = tomlInfo(resource, "names")
        print(text, file = tomlFile)

      # Write out the ops
      ops = "ops=["
      for i in range(0, noValues - 1): ops += "\"self\", "
      ops += "\"self\"]"
      print(ops, file = tomlFile)
      print(file = tomlFile)

# Include information in the toml file
def tomlInfo(resource, infoType):
  global resourceInfo
  noValues = len(resourceInfo["resources"][resource][infoType])

  # Define the fields to use
  text = infoType + "=["
  for index, value in enumerate(resourceInfo["resources"][resource][infoType]):
    if infoType == "columns": text += str(value)
    else: text += "\"" + str(value) + "\""
    if (index + 1) < noValues: text += ", "
  text += "]"

  return text, noValues

# Generate the command line for the annoatation script
def genBashScript(args):
  global resourceInfo
  global tsvFiles

  # Create a script file
  try: bashFile = open("calypso_annotation_pipeline.sh", "w")
  except: fail("There was a problem opening a file (calypso_annotation_pipeline.sh) to write to")

  # Initial information
  print("#! /bin/bash", file = bashFile)
  print("set -eou pipefail", file = bashFile)
  print(file = bashFile)

  # Define the names of the input and output files
  print("# Following are the input VCF and output files created by the pipeline", file = bashFile)
  vcfFile = os.path.abspath(args.input_vcf)
  print("VCF=", vcfFile, sep = "", file = bashFile)

  # Generate the names of the intermediate and final vcf files
  vcfBase = os.getcwd() + "/" + os.path.abspath(args.input_vcf).split("/")[-1].rstrip("vcf.gz")
  print("NORMVCF=" + str(vcfBase) + "_norm.vcf.gz", sep = "", file = bashFile)
  print("CLEANVCF=" + str(vcfBase) + "_clean.vcf.gz", sep = "", file = bashFile)
  print("ANNOTATEDVCF=" + str(vcfBase) + "_annotated.vcf.gz", sep = "", file = bashFile)
  if isFamily: print("SLIVAR1VCF=" + str(vcfBase) + "_slivar1.vcf.gz", sep = "", file = bashFile)
  if isFamily: print("SLIVAR2VCF=" + str(vcfBase) + "_slivar2.vcf.gz", sep = "", file = bashFile)
  finalVcf    = str(vcfBase) + "_calypso.vcf.gz"
  filteredVcf = str(vcfBase) + "_calypso_filtered.vcf.gz"
  print("FINALVCF=" + finalVcf, sep = "", file = bashFile)
  print("FILTEREDVCF=" + filteredVcf, sep = "", file = bashFile)

  # Write the ped file, if necessary
  if isFamily: print("PED=", os.path.abspath(args.ped), sep = "", file = bashFile)
  print(file = bashFile)

  # Define the required resources
  print("# Following is a list of required resources", file = bashFile)
  try: print("REF=", resourceInfo["resources"]["fasta"]["file"], sep = "", file = bashFile)
  except: fail("The resources json does not define a reference fasta file")
  try: print("GFF=", resourceInfo["resources"]["gff"]["file"], sep = "", file = bashFile)
  except: fail("The resources json does not define a gff file")
  try: print("GNOMAD=", resourceInfo["resources"]["gnomAD"]["file"], sep = "", file = bashFile)
  except: fail("The resources json does not define a gnomAD zip file")
  try: print("JS=", resourceInfo["resources"]["slivar_js"]["file"], sep = "", file = bashFile)
  except: fail("The resources json does not define the Slivar functions js file")
  print("TOML=", os.getcwd(), "/calypso_annotations.toml", sep = "", file = bashFile)
  print(file = bashFile)

  # Print out status messages
  print("# Normalize and subset original VCF", file = bashFile)
  print("  echo \"Starting annotation...\"", file = bashFile)
  print("  echo \"Subsetting and normalizing input VCF...\"", file = bashFile)
  print(file = bashFile)

  # Generate a samples text file from the ped file, if this is a family
  if isFamily: print("  tail -n+2 $PED | cut -f 2 | sort -u > samples.txt", file = bashFile)
  else: fail("Can't handle singletons yet. Need to build samples file in case there are background samples")
  print("  bcftools norm -m - -w 10000 -f $REF $VCF \\", file = bashFile)
  print("  | bcftools view -a -c 1 -S samples.txt -O z -o $NORMVCF", file = bashFile)
  print(file = bashFile)

  # Index the normalized vcf file
  print("# Index the normalized VCF file", file = bashFile)
  print("  bcftools index -t $NORMVCF", file = bashFile)
  print(file = bashFile)

  # Strip all existing annotations from the VCF file
  print("# Strip existing VCF annotations", file = bashFile)
  print("  echo \"Removing any existing annotations from input VCF...\"", file = bashFile)
  print(file = bashFile)
  print("  bcftools annotate -x INFO $NORMVCF -O z -o $CLEANVCF", file = bashFile)
  print("  bcftools index -t $CLEANVCF", file = bashFile)
  print(file = bashFile)

  # Delete the normalized vcf
  print("# Remove the normalized VCF file", file = bashFile)
  print("  rm -f $NORMVCF", file = bashFile)
  print("  rm -f \"$NORMVCF\".tbi", file = bashFile)
  print("  rm -f samples.txt", file = bashFile)
  print(file = bashFile)

  # Annotate the vcf
  print("# Annotate with bcftools csq and add further annotations with vcfanno", file = bashFile)
  print("  echo \"Annotating cleaned VCF...\"", file = bashFile)
  print(file = bashFile)
  print("  bcftools csq -f $REF --ncsq 40 -l -g $GFF $CLEANVCF \\", file = bashFile)
  print("  | vcfanno -p 16 $TOML /dev/stdin \\", file = bashFile)
  print("  | bcftools view -O z -o $ANNOTATEDVCF -", file = bashFile)
  print(file = bashFile)

  # Delete the clean vcf
  print("# Remove the clean VCF and toml files", file = bashFile)
  print("  rm -f $CLEANVCF", file = bashFile)
  print("  rm -f \"$CLEANVCF\".tbi", file = bashFile)
  print("  rm -f $TOML", file = bashFile)
  print(file = bashFile)

############### SINGLETONS NEED GNOMAD ANNOTATIONS FROM SLIVAR
  # If this is a family, run Slivar
  print("# Run Slivar", file = bashFile)
  print("  slivar_static expr \\", file = bashFile)
  print("  --vcf $ANNOTATEDVCF \\", file = bashFile)
  print("  --ped $PED \\", file = bashFile)
  print("  --js $JS \\", file = bashFile)
  print("  -g $GNOMAD \\", file = bashFile)
  print("  --info 'INFO.gnomad_popmax_af < 0.01 && variant.FILTER == \"PASS\" && variant.ALT[0] != \"*\"' \\", file = bashFile)
  print("  --family-expr 'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < 0.001' \\", file = bashFile)
  print("  --family-expr 'x_denovo:(variant.CHROM == \"X\" || variant.CHROM == \"chrX\") && fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af < 0.001' \\", file = bashFile)
  print("  --family-expr 'recessive:fam.every(segregating_recessive)' \\", file = bashFile)
  print("  --family-expr 'dominant:fam.every(segregating_dominant)' \\", file = bashFile)
  print("  -o $SLIVAR1VCF", file = bashFile)
  print("  bcftools index -t $SLIVAR1VCF", file = bashFile)
  print(file = bashFile)

  # Slivar compound hets
  print("# Use Slivar to determine compound hets", file = bashFile)
  print("  slivar_static expr \\", file = bashFile)
  print("  --pass-only \\", file = bashFile)
  print("  --vcf $ANNOTATEDVCF \\", file = bashFile)
  print("  --ped $PED \\", file = bashFile)
  print("  --js $JS \\", file = bashFile)
  print("  -g $GNOMAD \\", file = bashFile)
  print("  --family-expr 'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < 0.001' \\", file = bashFile)
  print("  --trio 'comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_popmax_af < 0.005' \\", file = bashFile)
  print("  | slivar_static compound-hets \\", file = bashFile)
  print("  --vcf /dev/stdin \\", file = bashFile)
  print("  --skip NONE \\", file = bashFile)
  print("  -s comphet_side \\", file = bashFile)
  print("  -s denovo \\", file = bashFile)
  print("  -p $PED \\", file = bashFile)
  print("  -o $SLIVAR2VCF", file = bashFile)
  print("  bcftools index -t $SLIVAR2VCF", file = bashFile)
  print(file = bashFile)

  # Delete the annotated vcf
  print("# Remove the annotated VCF", file = bashFile)
  print("  rm -f $ANNOTATEDVCF", file = bashFile)
  print("  rm -f \"$ANNOTATEDVCF\".tbi", file = bashFile)
  print(file = bashFile)

  # Combine the slivar vcf files
  print("# Combine the slivar vcf files", file = bashFile)
  print("  bcftools concat -a -d none $SLIVAR1VCF $SLIVAR2VCF -O z -o $FINALVCF", file = bashFile)
  print("  bcftools index -t $FINALVCF", file = bashFile)

  # Delete the slivar vcf files
  print("# Remove the Slivar VCF files", file = bashFile)
  print("  rm -f $SLIVAR1VCF", file = bashFile)
  print("  rm -f \"$SLIVAR1VCF\".tbi", file = bashFile)
  print("  rm -f $SLIVAR2VCF", file = bashFile)
  print("  rm -f \"$SLIVAR2VCF\".tbi", file = bashFile)
  print(file = bashFile)
  print("  echo \"Everything completed! Annotated VCF written to $FINALVCF\"", file = bashFile)
  print(file = bashFile)

  # Filter the VCF file to generate variants to pass to Mosaic
  print("# Filter the VCF file to generate variants to pass to Mosaic", file = bashFile)
  print("  echo \"Filtering final VCF file\"", file = bashFile)
  print("  slivar_static expr --vcf $FINALVCF \\", file = bashFile)
  print("  --info 'INFO.gnomad_popmax_af < 0.01 && variant.FILTER == \"PASS\" && variant.ALT[0] != \"*\"' \\", file = bashFile)
  print("  --pass-only \\", file = bashFile)
  print("  -o $FILTEREDVCF", file = bashFile)
  print(file = bashFile)

  # Generate the tsv files to pass annotations to Mosaic
  print("  # Generate the tsv files to pass annotations to Mosaic", file = bashFile)
  print("  echo \"Generating annotations tsv files for Mosaic\"", file = bashFile)
  print(file = bashFile)
  generateTsv(args, bashFile)

  # Close the file
  bashFile.close()

  # Make the annotation script executable
  makeExecutable = os.popen("chmod +x calypso_annotation_pipeline.sh").read()

  # Return the output vcf
  return finalVcf

# Generate the tsv files to pass annotations to Mosaic
def generateTsv(args, bashFile):
  global resourceInfo
  global mosaicInfo
  privateAnnotations = []

  # Loop over all the annotations to pass to Mosaic
  for resource in resourceInfo["resources"]:

    # Check if this resource is to be uploaded to Mosaic, and if it can be (e.g. does it appear
    # in the Mosaic json
    upload    = resourceInfo["resources"][resource]["upload"]
    canUpload = True if resource in mosaicInfo["resources"] else False

    # If the resource is listed as to be uploaded to Mosaic, but there are no instructions on how
    # to, fail
    if upload and not canUpload: fail("Resource '" + str(resource) + "' is marked as to be uploaded to Mosaic, but it is not included in the Mosaic json")

    # If the resource is not listed as to be uploaded to Mosaic, but it can be, provide a warning.
    elif not upload and canUpload:
      warningTitle = "Resource omitted"
      description  = "The following resources were listed as to be uploaded to Mosaic in the resources json, "
      description += "but no instructions are provided in the Mosaic json, so they will not be uploaded"
      if warningTitle not in warnings: warnings[warningTitle] = {"description": description, "messages": [resource]}
      else: warnings[warningTitle]["messages"].append(resource)

    # For resources to upload, check if they are public or private Mosaic annotations. All private
    # annotations should be included in the same tsv, but public annotation will get their own tsv.
    elif upload and canUpload:
      annotation_type = mosaicInfo["resources"][resource]["annotation_type"]
      if annotation_type == "private": privateAnnotations.append(resource)
      else: buildTsv(args, [resource], True, bashFile)

  # Build the tsv for private annotations
  buildTsv(args, privateAnnotations, False, bashFile)

# Build the tsv file
def buildTsv(args, resources, isPublic, bashFile):
  global mosaicInfo
  global tsvFiles

  # First create the header line
  header = "CHROM\\tSTART\\tEND\\tREF\\tALT"
  for resource in resources:
    for annotation in sorted(mosaicInfo["resources"][resource]["annotations"]):
      header += "\\t" + mosaicInfo["resources"][resource]["annotations"][annotation]["uid"]

  # Print the header
  if isPublic:
    print("  # Public annotation: ", resource, sep = "", file = bashFile)
    outputFile = resource + "_annotations_mosaic.tsv"
    tsvFiles[resource] = outputFile
  else:
    print("  # Private annotations", sep = "", file = bashFile)
    outputFile = "private_annotations_mosaic.tsv"
    tsvFiles["private"] = outputFile
  tempFile   = outputFile + ".tmp"
  print("  echo -e \"", header, "\" > ", outputFile, sep = "", file = bashFile)

  # Include all annotations for all resources
  annotations = "%CHROM\\t%POS\\t%END\\t%REF\\t%ALT"
  for resource in resources:
    for annotation in sorted(mosaicInfo["resources"][resource]["annotations"]):
      annotations += "\\t%INFO/" + annotation
  print("  bcftools query -f '", annotations, "\\n' \\", sep = "", file = bashFile)
  print("  $FILTEREDVCF >> ", outputFile, sep = "", file = bashFile)
  print(file = bashFile)
  print("  # Update the tsv file to Mosaic specifications", file = bashFile)
  print("  awk '{FS=\"\\t\"; OFS=\"\\t\"} {if ($1 != \"CHROM\") {if ($1 ~ \"chr\") {$1=substr($1, 4)}; $3 = $3+1; for(i=6; i<=NF; i++) {$i = ($i == \".\" ? \"\" : $i)}}; print $0}' \\", file = bashFile)
  print(" ", outputFile, ">", tempFile, file = bashFile)
  print("  mv", tempFile, outputFile, file = bashFile)
  print(file = bashFile)

# Parse the config file and get the Mosaic token and url
def parseConfig(args):
  global mosaicToken
  global mosaicUrl

  # Check the config file exists, if it was defined
  if args.config:
    if not exists(args.config): fail("Config file '" + str(args.config) + "' does not exist")

    # Parse the config file and store the token and url
    try: configFile = open(args.config, "r")
    except: fail("Failed to open config file '" + str(args.config) + "'")
    for line in configFile.readlines():
      fields = line.rstrip().split("=")
      if fields[0].startswith("MOSAIC_TOKEN"): mosaicToken = fields[1]
      elif fields[0].startswith("MOSAIC_URL"): mosaicUrl = fields[1]

    # Fail if the config file did not contain the required fields
    if not mosaicToken: fail("Config file '" + str(args.config) + "' does not contain the token (MOSAIC_TOKEN=***)")
    if not mosaicUrl: fail("Config file '" + str(args.config) + "' does not contain the url (MOSAIC_URL=***)")

    # Ensure the url terminates with a '/'
    if not mosaicUrl.endswith("/"): mosaicUrl += "/"

# Output a script to upload variants to Mosaic
def uploadVariants(args):
  global mosaicInfo

  # Open a script file
  uploadFileName  = "calypso_upload_variants_to_mosaic.sh"
  try: uploadFile = open(uploadFileName, "w")
  except: fail("Could not open " + str(uploadFileName) + " to write to")

  # Define the name of the variant set that will be create
  description = "Calypso_v" + str(mosaicInfo["version"]) + "_variants_" + str(date.today())

  # Write the command to file
  print("# Upload variants to Mosaic", file = uploadFile)
  print("python ", os.path.dirname( __file__ ), "/mosaic_commands/upload_annotations_to_mosaic.py \\", sep = "", file = uploadFile)
  print("  -i ", str(os.path.abspath(args.input_vcf)), " \\", sep = "", file = uploadFile)
  print("  -a ", os.path.dirname( __file__), "/mosaic_commands \\", sep = "", file = uploadFile)
  print("  -d \"", description, "\" \\", sep = "", file = uploadFile)
  if args.config: print("  -c", str(args.config), "\\", file = uploadFile)
  else: print("  -c \"Insert config file here\" \\", sep = "", file = uploadFile)
  if args.project_id: print("  -p", str(args.project_id), file = uploadFile)
  else: print("  -p \"Insert Mosaic project id here\"", file = uploadFile)

  # Close the file
  uploadFile.close()

  # Make the annotation script executable
  makeExecutable = os.popen("chmod +x " + str(uploadFileName)).read()

# Output scripts to upload annotations to Mosaic
def uploadAnnotations(args):
  global tsvFiles
  global mosaicInfo

  for resource in tsvFiles:

    # Open a script file
    uploadFileName  = "calypso_upload_" + str(resource) + "_annotations_to_mosaic.sh"
    try: uploadFile = open(uploadFileName, "w")
    except: fail("Could not open " + str(uploadFileName) + " to write to")

    # Write the command to file
    print("# Upload ", resource, " annotations to Mosaic", file = uploadFile)
    print("python ", os.path.dirname( __file__ ), "/mosaic_commands/upload_annotations_to_mosaic.py \\", sep = "", file = uploadFile)
    print("  -i ", str(tsvFiles[resource]), " \\", sep = "", file = uploadFile)
    print("  -a ", os.path.dirname( __file__), "/mosaic_commands \\", sep = "", file = uploadFile)
    if args.config: print("  -c", str(args.config), "\\", file = uploadFile)
    else: print("  -c \"Insert config file here\" \\", sep = "", file = uploadFile)

    # Public annotations need to be uploaded to the project that contains them, private
    # annotations are uploaded to the defined project.
    if resource == "private":
      if args.project_id: print("  -p", str(args.project_id), file = uploadFile)
      else: print("  -p \"Insert Mosaic project id here\"", file = uploadFile)
    else: print("  -p ", mosaicInfo["resources"][resource]["project_id"], sep = "", file = uploadFile)

    # Close the file
    uploadFile.close()

    # Make the annotation script executable
    makeExecutable = os.popen("chmod +x " + str(uploadFileName)).read()

# Output summary file
def calypsoSummary(args, finalVcf):
  global resourceInfo
  global version
  today = date.today()

  # Open an output summary file
  try: summaryFile = open("calypso_" + str(today) + ".txt", "w")
  except: fail("Failed to open summary file")

  # Write relevant information to file
  print("### Output from Calypso pipeline ###", file = summaryFile)
  print(file = summaryFile)
  print("Calypso pipeline version: ", version, sep = "", file = summaryFile)
  print("Created on:               ", today, sep = "", file = summaryFile)
  print("Generated VCF file:       ", finalVcf, sep = "", file = summaryFile)
  print(file = summaryFile)

  # Loop over all the used resources and output their versions
  for resource in resourceInfo["resources"]:
    print(resource, file = summaryFile)
    print("  version: ", resourceInfo["resources"][resource]["version"], sep = "", file = summaryFile)
    print("  file:    ", resourceInfo["resources"][resource]["file"], sep = "", file = summaryFile)
    print(file = summaryFile)

  # Close the file
  summaryFile.close()

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

# Pipeline version
version = "0.1.1"

# Store warnings to be output at the end
warnings = {}

# Store info on allowed values
isFamily          = True
allowedFamily     = ['singleton', 'duo', 'trio']
allowedReferences = ['37', '38']

# Store information about the resource files
resourceVersion = False
resourceInfo    = {}

# Store information related to Mosaic
mosaicToken = False
mosaicUrl   = False
mosaicInfo  = {}

# Store the tsv files that upload annotations to Mosaic
tsvFiles = {}

if __name__ == "__main__":
  main()
