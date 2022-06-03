#!/usr/bin/python

from __future__ import print_function
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

  # Build the toml file for vcfanno
  buildToml()

  # Generate the bash script to run the annotation pipeline
  genBashScript(args)

# Input options
def parseCommandLine():
  parser = argparse.ArgumentParser(description='Process the command line')
  parser.add_argument('--reference', '-r', required = True, metavar = "string", help = "The reference genome to use. Allowed values: '37', '38'")
  parser.add_argument('--family_type', '-f', required = True, metavar = "string", help = "The familty structure. Allowed values: 'singleton', 'duo', 'trio'")
  parser.add_argument('--data', '-d', required = True, metavar = "string", help = "The path to the directory where the resources live")
  parser.add_argument('--vcf', '-i', required = True, metavar = "string", help = "The input vcf file to annotate")
  parser.add_argument('--ped', '-p', required = False, metavar = "string", help = "The pedigree file for the family. Not required for singletons")
  parser.add_argument('--version', '-v', action="version", version='%(prog)s 0.1.1')

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
  if not exists(args.vcf): fail("The vcf file could not be found")
  elif not args.vcf.endswith(".vcf.gz"): fail("The input vcf file (--vcf, -v) must be a compressed, indexed vcf and have the extension '.vcf.gz'")
  if isFamily:
    if not args.ped: fail("A ped file needs to specified (--ped, -p) for family type \"" + str(args.family_type) + "\"")
    if not exists(args.ped): fail("The ped file could not be found")
    elif not args.ped.endswith(".ped"): fail("The input ped file (--ped, -p) must have the extension '.ped'")

# Parse the config file describing the resources for the selected genome build,
# check the files exist, and store the versions.
def checkResources(args):
  global resourceInfo

  # Ensure the data path ends with a "/", then add the refrence directory
  if args.data[-1] != "/": args.data += "/"
  resourceInfo["path"] = args.data + "GRCh" + str(args.reference) + "/"

  # Define the name of the resource description json file
  if args.reference == '37': resourceFilename = args.data + "resources_GRCh37.json"
  elif args.reference == '38': resourceFilename = args.data + "resources_GRCh38.json"

  # Try and open the file
  try: resourceFile = open(resourceFilename, "r")
  except: fail("The file describing the resource files (" + str(resourceFilename) + ") could not be found in the data directory")

  # Extract the json information
  try: resourceData = json.loads(resourceFile.read())
  except: fail("Not valid json")

  # Store the data version
  try: resourceInfo["version"] = resourceData['version']
  except: fail("The resource json does not include a version")

  # Check that the resource json reference matches the selected reference
  try: resourceReference = resourceData['reference']
  except: fail("The resource json does not include a reference genome")
  isRefMatch = False
  if args.reference == '37' and resourceReference == 'GRCh37': isRefMatch = True
  elif args.reference == '38' and resourceReference == 'GRCh38': isRefMatch = True
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
    except: fail("The resources json did not indicate if resource \"" + str(resource) + "\" should be included in the toml file")

    # If the resource is to be included in a toml file, the fields in the vcf INFO need to be specified
    if resourceInfo["resources"][resource]["toml"]:
      if "fields" in resources[resource]: resourceInfo["resources"][resource]["fields"] = resources[resource]["fields"]
      if "columns" in resources[resource]: resourceInfo["resources"][resource]["columns"] = resources[resource]["columns"]
      if "names" in resources[resource]: resourceInfo["resources"][resource]["names"] = resources[resource]["names"]

    # Check that the file exists
    if not exists(resourceInfo["resources"][resource]["file"]): fail("Resource file " + str(resourceInfo["resources"][resource]["file"]) + " does not exist")

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

  # Create a script file
  try: bashFile = open("calypso_annotation_pipeline.sh", "w")
  except: fail("There was a problem opening a file (calypso_annotation_pipeline.sh) to write to")

  # Initial information
  print("#! /bin/bash", file = bashFile)
  print("set -eou pipefail", file = bashFile)
  print(file = bashFile)

  # Define the names of the input and output files
  print("# Following are the input VCF and output files created by the pipeline", file = bashFile)
  vcfFile = os.path.abspath(args.vcf)
  print("VCF=", vcfFile, sep = "", file = bashFile)

  # Generate the names of the intermediate and final vcf files
  vcfBase = os.getcwd() + "/" + os.path.abspath(args.vcf).split("/")[-1].rstrip("vcf.gz")
  print("NORMVCF=" + str(vcfBase) + "_norm.vcf.gz", sep = "", file = bashFile)
  print("CLEANVCF=" + str(vcfBase) + "_clean.vcf.gz", sep = "", file = bashFile)
  print("ANNOTATEDVCF=" + str(vcfBase) + "_annotated.vcf.gz", sep = "", file = bashFile)
  if isFamily: print("SLIVAR1VCF=" + str(vcfBase) + "_slivar1.vcf.gz", sep = "", file = bashFile)
  if isFamily: print("SLIVAR2VCF=" + str(vcfBase) + "_slivar2.vcf.gz", sep = "", file = bashFile)
  print("FINALVCF=" + str(vcfBase) + "_calypso.vcf.gz", sep = "", file = bashFile)

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

  # Make the annotation script executable
  makeExecutable = os.popen("chmod +x calypso_annotation_pipeline.sh").read()

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)

# Initialise global variables
isFamily          = True
allowedFamily     = ['singleton', 'duo', 'trio']
allowedReferences = ['37', '38']

# Store information about the resource files
resourceVersion = False
resourceInfo    = {}

if __name__ == "__main__":
  main()
