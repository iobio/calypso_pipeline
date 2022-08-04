#!/usr/bin/python

from __future__ import print_function
from datetime import date
from os.path import exists

import os
import argparse
import json
import math
import glob

# Add the path of the common functions and import them
from sys import path
path.append(os.path.dirname(__file__) + "/mosaic_commands")
import mosaic_config
import api_variant_annotations as api_va
#import api_projects as api_p
import api_project_attributes as api_pa
import api_samples as api_s
import api_variant_filters as api_vf

def main():
  global mosaicConfig

  # Parse the command line
  args = parseCommandLine()

  # Check the supplied parameters are as expected
  checkArguments(args)

  # Check all resources
  checkResources(args)

  # Create a directory to store all created files
  setWorkingDir()

  # Parse the Mosaic config file to get the token and url for the api calls
  mosaicRequired = {"token": True, "url": True, "attributeProjectId": False}
  mosaicConfig   = mosaic_config.parseConfig(args, mosaicRequired)

  # Parse the Mosaic json
  parseMosaicJson(args)

  # If a previous Calypso vcf file was set, get the version assigned to the vcf
  if args.previous_vcf: getPreviousVcfInfo(args)

  # Determine the id of the proband
  if args.ped: getProband(args)

  # Get the order that the samples appear in, in the vcf header
  getSampleOrder(args)

  # Build the toml file for vcfanno
  buildToml()

  # Generate the bash script to run the annotation pipeline
  finalVcf, filteredVcf = genBashScript(args)

  # Generate script to upload filtered variants to Mosaic
  uploadVariants(args, filteredVcf)

  # Generate scripts to upload annotations
  uploadAnnotations(args)

  # Create Mosaic filters
  variantFilters(args)

  # Output summary file
  calypsoSummary(args, finalVcf)

  # Set project attributes to indicate when and which version of Calypso has been run
  updateCalypsoAttributes(args)

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
  parser.add_argument('--project_id', '-p', required = True, metavar = "string", help = "The project id that variants will be uploaded to")

  # Optional reanalysis arguments
  parser.add_argument('--previous_vcf', '-s', required = False, metavar = "string", help = "A previously annotated Calypso vcf file")

  # Optional pipeline arguments
  parser.add_argument('--ped', '-e', required = False, metavar = "string", help = "The pedigree file for the family. Not required for singletons")
  parser.add_argument('--resource_json', '-j', required = False, metavar = "string", help = "The json file describing the annotation resources")

  # Optional mosaic arguments
  parser.add_argument('--config', '-c', required = False, metavar = "string", help = "The config file for Mosaic")
  parser.add_argument('--token', '-t', required = False, metavar = "string", help = "The Mosaic authorization token")
  parser.add_argument('--url', '-u', required = False, metavar = "string", help = "The base url for Mosaic")
  parser.add_argument('--attributes_project', '-a', required = False, metavar = "integer", help = "The Mosaic project id that contains public attributes")
  parser.add_argument('--mosaic_json', '-m', required = False, metavar = "string", help = "The json file describing the Mosaic parameters")

  # Version
  parser.add_argument('--version', '-v', action="version", version='Calypso annotation pipeline version: ' + str(version))

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
  global workingDir

  # Ensure the data path ends with a "/", then add the refrence directory
  if args.data_directory[-1] != "/": args.data_directory += "/"
  resourceInfo["path"] = args.data_directory + "GRCh" + str(args.reference) + "/"

  # Define the name of the resource description json file. Use the file provided on the command
  # line, and resort to the default file if not included. The version of the resource file is
  # unknown, so look for any json with the correct prefix and if there is more than one such
  # file, throw an error.
  if args.resource_json: resourceFilename = args.resource_json
  else:
    if args.reference == '37': resourceFilename = args.data_directory + "resources_GRCh37"
    elif args.reference == '38': resourceFilename = args.data_directory + "resources_GRCh38"
    resourceFiles = glob.glob(resourceFilename + "*json")
    if len(resourceFiles) != 1: fail("There are zero, or more than one resource files for GRCh" + args.reference + " in " + args.data_directory)
    resourceFilename = resourceFiles[0]

  # Try and open the file
  try: resourceFile = open(resourceFilename, "r")
  except: fail("The file describing the resource files (" + str(resourceFilename) + ") could not be found")

  # Extract the json information
  try: resourceData = json.loads(resourceFile.read())
  except: fail("The json file (" + str(resourceFilename) + ") is not valid")

  # Store the resource data version
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
    try:
      resourceInfo["resources"][resource]["file"] = resources[resource]["file"]
      if resourceInfo["resources"][resource]["file"]: resourceInfo["resources"][resource]["file"] = resourceInfo["path"] + resources[resource]["file"]
    except: fail("File for resource \"" + str(resource) + "\" was not included in the resources json")
    try: resourceInfo["resources"][resource]["toml"] = resources[resource]["toml"]
    except: fail("The resources json did not indicate if resource \"" + str(resource) + "\" should be included in the toml file")
    try: resourceInfo["resources"][resource]["upload"] = resources[resource]["upload_to_mosaic"]
    except: fail("The resources json did not indicate if resource \"" + str(resource) + "\" should be uploaded to Mosaic")

    # Only check for the "ops" field if this is to be included in the toml file
    if resourceInfo["resources"][resource]["toml"]:
      try: resourceInfo["resources"][resource]["ops"] = resources[resource]["ops"]
      except: fail("The resources json did not indicate the toml \"ops\" commands for resource \"" + str(resource) + "\"")

    # If the resource is to be included in a toml file, the fields in the vcf INFO need to be specified
    if resourceInfo["resources"][resource]["toml"]:
      if "fields" in resources[resource]: resourceInfo["resources"][resource]["fields"] = resources[resource]["fields"]
      if "columns" in resources[resource]: resourceInfo["resources"][resource]["columns"] = resources[resource]["columns"]
      if "names" in resources[resource]: resourceInfo["resources"][resource]["names"] = resources[resource]["names"]

    # Check that the file exists
    if not exists(resourceInfo["resources"][resource]["file"]): fail("Resource file " + str(resourceInfo["resources"][resource]["file"]) + " does not exist")

    # If the resource is vep, handle this separately
    if resource == "vep": processVep(resources)

    # If this is a resource that uses VEP to generate, process the VEP fields
    if "vep_commands" in resources[resource]: processVepCommands(resources, resource)
    else: resourceInfo["resources"][resource]["isVepAnnotation"] = False

    # Check if this resource is to be applied to the full file, or only applied to the filtered file. If the
    # annotation takes a long time to process, and does not need to be tracked longitudinally (e.g. HGVS), the
    # annotation can be applied to the filtered file only, to generate values to upload to Mosaic.
    isFiiltered = resources[resource]["apply_to_filtered"] = True if "apply_to_filtered" in resources[resource] else False
    resourceInfo["resources"][resource]["apply_to_filtered"] = isFiiltered

    # Check if there are any postannotation commands for vcfanno
    if "post_annotation" in resources[resource]:
      resourceInfo["resources"][resource]["post_annotation"] = {}
      for postAnnoField in resources[resource]["post_annotation"]:
        if postAnnoField == "name": resourceInfo["resources"][resource]["post_annotation"]["name"] = resources[resource]["post_annotation"]["name"]
        elif postAnnoField == "fields": resourceInfo["resources"][resource]["post_annotation"]["fields"] = resources[resource]["post_annotation"]["fields"]
        elif postAnnoField == "op": resourceInfo["resources"][resource]["post_annotation"]["op"] = resources[resource]["post_annotation"]["op"]
        elif postAnnoField == "type": resourceInfo["resources"][resource]["post_annotation"]["type"] = resources[resource]["post_annotation"]["type"]
        else: fail("Unexpected post_annotation field in the resources json for resource \"" + str(resource) + "\"")

# Create a working directory
def setWorkingDir():
  global resourceInfo
  global workingDir

  workingDir += resourceInfo['version'] + "/"

  # If the directory doesn't exist, create it
  if not os.path.exists(workingDir): os.makedirs(workingDir)
  
# Process the vep information
def processVep(resources):
  global resourceInfo

  # VEP requires a cache and plugins directory to run. Get these directories and check they exist
  try: resourceInfo["resources"]["vep"]["cache"] = resources["vep"]["cache"]
  except: fail("The VEP resource description does not include the \"cache\"")
  try: resourceInfo["resources"]["vep"]["plugins"] = resources["vep"]["plugins"]
  except: fail("The VEP resource description does not include \"plugins\"")

  if not exists(resourceInfo["resources"]["vep"]["cache"]): fail("VEP cache does not exist")
  if not exists(resourceInfo["resources"]["vep"]["plugins"]): fail("VEP plugins directory does not exist")

# For resources that use VEP to generate annotations, process the VEP information
def processVepCommands(resources, resource):
  global resourceInfo

  # Get the VEP specific information for this resource
  resourceInfo["resources"][resource]["isVepAnnotation"] = True
  if "commands" in resources[resource]["vep_commands"]: resourceInfo["resources"][resource]["vep_commands"] = resources[resource]["vep_commands"]["commands"]
  if "fields" in resources[resource]["vep_commands"]: resourceInfo["resources"][resource]["vep_fields"] = resources[resource]["vep_commands"]["fields"]

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

    # Some annotations will require postprocessing with specific tools. Store this information
    mosaicInfo["resources"][resource]["postprocess"] = False
    if "postprocess" in resources[resource]:
      mosaicInfo["resources"][resource]["postprocess"] = resources[resource]["postprocess"]
      if "postprocess_command" not in resources[resource]: fail("Resource \"" + resource + "\" requires postprocessing, but no tool has been specified")

      # Store information on the postprocessing command
      try: precommand = resources[resource]["postprocess_command"]["precommand"]
      except: fail("Resource \"" + resource + "\" postprocess_command is missing \"precommand\"")
      try: tool = resources[resource]["postprocess_command"]["tool"]
      except: fail("Resource \"" + resource + "\" postprocess_command is missing \"tool\"")
      try: splitTie = resources[resource]["postprocess_command"]["split_tie"]
      except: fail("Resource \"" + resource + "\" postprocess_command is missing \"split_tie\"")
      try: arguments = resources[resource]["postprocess_command"]["args"]
      except: fail("Resource \"" + resource + "\" postprocess_command is missing \"args\"")
      mosaicInfo["resources"][resource]["post_command"] = {"precommand": precommand, "tool": tool, "args": arguments, "split_tie": splitTie}

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

      # Some annotations are pulled directly from native VCF genotype fields (e.g. FORMAT). This will be indicated with
      # the "isGenotype" field.
      isGenotype = resources[resource]["annotations"][annotation]["is_genotype"] if "is_genotype" in resources[resource]["annotations"][annotation] else False

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
            mosaicInfo["resources"][resource]["annotations"][annotation] = {"uid": uid, "type": annType, "operation": operation, "positions": maxPositions, "isGenotype": isGenotype}
          else: fail("Annotation '" + str(annotation) + "' for '" + str(resource) + "' include \"positions\" defining the fields to find the max of, but this needs to be a list")

      # If the delimiter is set, a "position" field must be present with a numerical value.
      else:
        if mosaicInfo["resources"][resource]["delimiter"]:
          try: position = resources[resource]["annotations"][annotation]["position"]
          except: fail("Annotation '" + str(annotation) + "' for '" + str(resource) + "' requires the \"position\" field as the delimiter is set")
          if not isinstance(position, int): fail("Annotation '" + annotation + "' for '" + str(resource) + "' has a position field that is not an integer")
          mosaicInfo["resources"][resource]["annotations"][annotation] = {"uid": uid, "type": annType, "position": position, "isGenotype": isGenotype}
        else: mosaicInfo["resources"][resource]["annotations"][annotation] = {"uid": uid, "type": annType, "isGenotype": isGenotype}

# If a previous Calypso vcf file was set, get the version assigned to the vcf
def getPreviousVcfInfo(args):
  global resourceInfo
  global workingDir
  global version
  headerVcf = workingDir + "calypso_previous_header.vcf"
  previousVersions = []
  previousFiles    = []

  # Check that the previous vcf file exists
  if not exists(args.previous_vcf): fail("Previous Calypso vcf file, " + args.previous_vcf + ", does not exist.")

  # Grab the vcf header and extract the previous versions
  try: os.popen("bcftools view -h -O v -o " + headerVcf + " " + args.previous_vcf).read()
  except: fail("Unable to extract header from vcf file, " + args.previous_vcf)
  try: headerInfo = open(headerVcf, "r")
  except: fail("Unable to read header information from file " + arg.previous_vcf)
  for line in  headerInfo.readlines():
    if line.rstrip().startswith("##calypsoVersion="): previousVersions.append(line.rstrip().split("=")[1])
    if line.rstrip().startswith("##calypsoSourceVcf="): previousFiles.append(line.rstrip().split("=")[1])

  # If there are multiple Calypso versions in the header, instruct the user to clean this up to only include the most recent
  if len(previousVersions) != 1: fail("The previous Calypso vcf file, " + args.previous_vcf + ", either has no, or multiple \"##calypsoVersion=VERSION\" header lines. Ensure it only has one")
  if len(previousFiles) != 1: fail("The previous Calypso vcf file, " + args.previous_vcf + ", either has no, or multiple \"##calypsoSourceVcf=VCFFILE\" header lines. Ensure it only has one")

  # Get the previous Calypso and the resource version
  versions = previousVersions[0].split("r")
  resourceInfo['previous_version'] = versions[1]
  previousVersion = versions[0].strip('v')
  previousFile    = previousFiles[0].split('/')[-1]

  # If both the Calypso and the resource versions are the same as for the previously generated vcf, and the source
  # vcf files are the same, no changes have been made, so there is no need to rerun the pipeline
  if str(resourceInfo['version']) == str(resourceInfo['previous_version']) and str(previousVersion) == str(version) and str(args.input_vcf) == str(previousFile):
    failString  = "Running the same Calypso version (" + version + "), with the same resource version (" + resourceInfo['version'] + "), on the same\nsource vcf ("
    failString += args.input_vcf + ") will result in the same output, so Calypso will not be run"
    fail(failString)

  # Delete the created header file
  try: os.remove(headerVcf)
  except: fail("Failed to delete created header vcf, " + headerVcf)

# Determine the id of the proband
def getProband(args):
  global mosaicConfig
  global samples
  global proband

  # Open the ped file and get the samples
  try: pedFile = open(args.ped, "r")
  except: fail("Couldn't open the ped file (" + args.ped + ")")

  # Get the samples Mosaic id and store
  mosaicSamples = {}
  for sample in json.loads(os.popen(api_s.getSamples(mosaicConfig, args.project_id)).read()):
    mosaicSamples[sample["name"]] = sample["id"]

  # Get information on all the samples
  noAffected = 0
  for line in pedFile:

    # Ignore the header line
    if not line.startswith("#"):
      fields   = line.rstrip().split("\t")
      sample   = fields[1]

      # Determine if the sample has parents
      father = False if fields[2] == "0" else fields[2]
      mother = False if fields[3] == "0" else fields[3]
      if mother and father: noParents = 2
      elif not mother and not father: noParents = 0
      else: noParents = 1

      # Determine if the sample is affected
      if int(fields[5]) == 2:
        isAffected = True 
        noAffected += 1
        proband = sample
      else: isAffected = False
      samples[sample] = {"noParents": noParents, "father": father, "mother": mother, "isAffected": isAffected, "relationship": False}
      if isAffected: samples[sample]["relationship"] = "Proband"

      # Attach the Mosaic sample id
      try: samples[sample]["mosaicId"] = mosaicSamples[sample]
      except: fail("Sample " + str(sample) + " is not present in the Mosaic project")

  # Check that the ped file conforms to the supplied family type
  if args.family_type == "singleton":
    if len(samples) != 1: fail("Family type was specified as a singleton, but the ped file contains multiple samples")
  elif args.family_type == "duo":
    if len(samples) != 2: fail("Family type was specified as a duo, but the ped file doesn't contain two samples")
  elif args.family_type == "trio":
    if len(samples) != 3: fail("Family type was specified as a trio, but the ped file doesn't contain three samples")
  elif args.family_type == "quad":
    if len(samples) != 4: fail("Family type was specified as a quad, but the ped file doesn't contain three samples")

  # If multiple samples are affected, throw a warning
  if noAffected != 1: fail("Cannot determine proband")

  # Identify the mother and father of the proband
  mother = samples[proband]["mother"]
  father = samples[proband]["father"]
  if samples[proband]["mother"]: samples[mother]["relationship"] = "Mother"
  if samples[proband]["father"]: samples[father]["relationship"] = "Father"

  # Identify siblings
  for sample in samples:
    if samples[sample]["mother"] == mother and samples[sample]["father"] and not samples[sample]["isAffected"]: samples[sample]["relationship"] = "Sibling"
    if not samples[sample]["relationship"]: fail("Sample " + str(sample) + " has an unknown relationship to the proband")

# Create variant filters in Mosaic
def variantFilters(args):
  global mosaicConfig
  global samples
  global rootPath
  global genotypeOptions
  global createAnnotations

  # Get information about the samples
  probandId = False
  motherId  = False
  fatherId  = False
  jsonPath  = rootPath + str("/mosaic_filters/")

  # Get the id of the proband and the parents
  for sample in samples:
    if samples[sample]["relationship"] == "Proband": probandId = samples[sample]["mosaicId"]
    elif samples[sample]["relationship"] == "Mother": motherId = samples[sample]["mosaicId"]
    elif samples[sample]["relationship"] == "Father": fatherId = samples[sample]["mosaicId"]

  # Get all of the filters that exist in the project
  existingFilters = {}
  try:
    for existingFilter in json.loads(os.popen(api_vf.getVariantFilters(mosaicConfig, args.project_id)).read()): existingFilters[existingFilter["name"]] = existingFilter["id"]
  except: fail("Unable to get existing variant filters for project " + str(args.project_id))
  
  # Loop over all the filter json files
  for filterJson in os.listdir(jsonPath):
    if filterJson.endswith(".json"):

      try: filterFile = open(jsonPath + filterJson, "r")
      except: fail("File, " + jsonPath + filterJson + ", could not be opened")
      try: data = json.load(filterFile)
      except: fail("File, " + jsonPath + filterJson + ", is not a valid json file")
      filterFile.close()

      # Get the name of the filter to create
      try: filterName = data["name"]
      except: fail("Varianr filter file, " + str(filterJson) + ", defines a filter with no name. A name needs to be supplied")

      # Check if the filter already exists in the project. If so, remove it
      if filterName in existingFilters:
        try: deleteFilter = os.popen(api_vf.deleteVariantFilter(mosaicConfig, args.project_id, existingFilters[filterName])).read()
        except: fail("Unable to delete variant filter, " + str(filterName) + " from project " + str(args.project_id))

      # Check what genotype filters need to be applied
      try: genotypes = data["genotypes"]
      except: fail("Mosaic variant filter json file, " + str(filterJson) + ", does not contain the \"genotypes\" field")
      for geno in genotypes:
        if geno not in genotypeOptions: fail("Mosaic variant filter json file, " + str(filterJson) + ", contains unknown genotype options: " + str(geno))
        if not genotypes[geno]: continue

        # Check which samples need to have the requested genotype and add to the command
        sampleList = []
        for sample in genotypes[geno]:
          if sample == "proband": sampleList.append(probandId)
          elif sample == "mother":
            if motherId: sampleList.append(motherId)
          elif sample == "father":
            if fatherId: sampleList.append(fatherId)
          elif sample == "sibling":
            fail("Can't handle siblings in filters yet")
          else: fail("Unknown sample, " + str(sample) + " in genotypes for Mosaic variant filter (" + str(filterJson) + ")")

        # Add the genotype filter to the filters listed in the json
        data["filters"][geno] = sampleList

      # Check the filters provided in the json. The annotation filters need to extract the
      # uids for the annotations
      removeAnnotations = []
      for index, annotationFilter in enumerate(data["filters"]["annotation_filters"]):

###################
###################
################### ADD CHECKS ON THE FILTER JSON
###################
###################

        # If the uid is provided with the filter, this annotation is complete. If not, the name
        # field must be present, and this must point to a created annotation uid. Remove the name
        # field and replace it with the uid.
        if "uid" not in annotationFilter:
          try: annotationName = annotationFilter["name"]
          except: fail("Annotation " + str(filterJson) + " does not have a uid, but also no name")
          annotationFilter.pop("name")

          # If the annotation is listed as optional and the name doesn't point to a created annotation
          # delete this annotation from the filter. This could be, for example, a filter on genotype
          # quality for the mother, but the mother isn't present in the project
          isOptional = annotationFilter.pop("optional") if "optional" in annotationFilter else False
          try: annotationFilter["uid"] = createdAnnotations[annotationName]
          except: 
            if isOptional: removeAnnotations.append(index)
            else: fail("Annotation, " + str(annotationName) + ", in " + str(filterJson) + " does not have a uid, and the name does not point to a created annotation")

      # Remove any optional filters that did not have uids
      for index in reversed(removeAnnotations): data["filters"]["annotation_filters"].pop(index)

      # Post a new filter
      execute = json.loads(os.popen(api_vf.postVariantFilter(mosaicConfig, data["name"], data["filters"], args.project_id)).read())

# Get the order that the samples appear in, in the vcf header
def getSampleOrder(args):
  global samples
  global sampleOrder

  command = "bcftools view -h " + str(args.input_vcf) + " | tail -1"
  data    = os.popen(command).read().rstrip().split("\t")
  for index, sample in enumerate(data[9:]):
    if sample in samples: sampleOrder.append(sample)

  # Check that every sample in samples is in sampleOrder
  for sample in samples:
    if sample not in sampleOrder: fail("Sample " + str(sample) + "is listed in the ped file, but does not appear in the vcf header")

# Build the toml file
def buildToml():
  global resourceInfo
  global workingDir
  global tomlFilename

  # Create a toml file
  tomlFilename  = workingDir + "calypso_annotations.toml"
  try: tomlFile = open(tomlFilename, "w")
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
      ops   = 'ops=["'
      noOps = len(resourceInfo["resources"][resource]["ops"])
      for i in range(0, noOps - 1): ops += str(resourceInfo["resources"][resource]["ops"][i]) + '", "'
      ops += str(resourceInfo["resources"][resource]["ops"][-1]) + '"]'
      print(ops, file = tomlFile)
      print(file = tomlFile)

      # If there is post annotation information to include in the toml, include it
      if "post_annotation" in resourceInfo["resources"][resource]:
        post = resourceInfo["resources"][resource]["post_annotation"]
        print("[[postannotation]]", file = tomlFile)
        if "name" in post: print('name="', post["name"], '"', sep = "", file = tomlFile)
        if "fields" in post:
          fieldString = 'fields=["'
          for i in range(0, len(post["fields"]) - 1): fieldString += str(post["fields"][i]) + '", "'
          fieldString += str(post["fields"][-1]) + '"]'
          print(fieldString, file = tomlFile)
        if "op" in post: print('op="', post["op"], '"', sep = "", file = tomlFile)
        if "type" in post: print('type="', post["type"], '"', sep = "", file = tomlFile)
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
  global workingDir
  global tsvFiles
  global tomlFilename

  # Create a script file
  bashFilename  = workingDir + "calypso_annotation_pipeline.sh"
  try: bashFile = open(bashFilename, "w")
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
  vcfBase = workingDir + os.path.abspath(args.input_vcf).split("/")[-1].rstrip("vcf.gz")
  print("NORMVCF=" + str(vcfBase) + "_norm.vcf.gz", sep = "", file = bashFile)
  print("CLEANVCF=" + str(vcfBase) + "_clean.vcf.gz", sep = "", file = bashFile)
  print("ANNOTATEDVCF=" + str(vcfBase) + "_annotated.vcf.gz", sep = "", file = bashFile)
  if isFamily: print("SLIVAR1VCF=" + str(vcfBase) + "_slivar1.vcf.gz", sep = "", file = bashFile)
  if isFamily: print("SLIVAR2VCF=" + str(vcfBase) + "_slivar2.vcf.gz", sep = "", file = bashFile)
  finalVcf    = str(vcfBase) + "_calypso.vcf.gz"
  filteredVcf = str(vcfBase) + "_calypso_filtered.vcf.gz"
  print("FINALVCF=" + finalVcf, sep = "", file = bashFile)
  print("FILTEREDVCF=" + filteredVcf, sep = "", file = bashFile)
  print("STDOUT=calypso_annotation_pipeline.stdout", file = bashFile)
  print("STDERR=calypso_annotation_pipeline.stderr", file = bashFile)

  # Write the ped file, if necessary
  if isFamily: print("PED=", os.path.abspath(args.ped), sep = "", file = bashFile)
  print(file = bashFile)

  # Define the required resources
  print("# Following is a list of required resources", file = bashFile)
  try: print("REF=", resourceInfo["resources"]["fasta"]["file"], sep = "", file = bashFile)
  except: fail("The resources json does not define a reference fasta file")
  try: print("GFF=", resourceInfo["resources"]["gff"]["file"], sep = "", file = bashFile)
  except: fail("The resources json does not define a gff file")
  try: print("SLIVAR_GNOMAD=", resourceInfo["resources"]["slivar_gnomAD"]["file"], sep = "", file = bashFile)
  except: fail("The resources json does not define a gnomAD zip file")
  try: print("JS=", resourceInfo["resources"]["slivar_js"]["file"], sep = "", file = bashFile)
  except: fail("The resources json does not define the Slivar functions js file")
  print("TOML=", tomlFilename, sep = "", file = bashFile)
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
  print("  2> $STDERR \\", file = bashFile)
  print("  | bcftools view -a -c 1 -S samples.txt -O z -o $NORMVCF \\", file = bashFile)
  print("  > $STDOUT 2>> $STDERR", file = bashFile)
  print(file = bashFile)

  # Index the normalized vcf file
  print("# Index the normalized VCF file", file = bashFile)
  print("  bcftools index -t $NORMVCF \\", file = bashFile)
  print(file = bashFile)

  # Strip all existing annotations from the VCF file
  print("# Strip existing VCF annotations", file = bashFile)
  print("  echo \"Removing any existing annotations from input VCF...\"", file = bashFile)
  print(file = bashFile)
  print("  bcftools annotate -x INFO $NORMVCF -O z -o $CLEANVCF \\", file = bashFile)
  print("  >> $STDOUT 2>> $STDERR", file = bashFile)
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
  print("  2>> $STDERR \\", file = bashFile)
  print("  | vcfanno -p 16 $TOML /dev/stdin \\", file = bashFile)
  print("  2>> $STDERR \\", file = bashFile)
  print("  | bcftools view -O z -o $ANNOTATEDVCF - \\", file = bashFile)
  print("  >> $STDOUT 2>> $STDERR", file = bashFile)
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
  print("  -g $SLIVAR_GNOMAD \\", file = bashFile)
  print("  --info 'INFO.gnomad_popmax_af < 0.01 && variant.FILTER == \"PASS\" && variant.ALT[0] != \"*\"' \\", file = bashFile)
  print("  --family-expr 'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < 0.001' \\", file = bashFile)
  print("  --family-expr 'x_denovo:(variant.CHROM == \"X\" || variant.CHROM == \"chrX\") && fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af < 0.001' \\", file = bashFile)
  print("  --family-expr 'recessive:fam.every(segregating_recessive)' \\", file = bashFile)
  print("  --family-expr 'dominant:fam.every(segregating_dominant)' \\", file = bashFile)
  print("  -o $SLIVAR1VCF \\", file = bashFile)
  print("  >> $STDOUT 2>> $STDERR", file = bashFile)
  print("  bcftools index -t $SLIVAR1VCF", file = bashFile)
  print(file = bashFile)

  # Slivar compound hets
  print("# Use Slivar to determine compound hets", file = bashFile)
  print("  slivar_static expr \\", file = bashFile)
  print("  --pass-only \\", file = bashFile)
  print("  --vcf $ANNOTATEDVCF \\", file = bashFile)
  print("  --ped $PED \\", file = bashFile)
  print("  --js $JS \\", file = bashFile)
  print("  -g $SLIVAR_GNOMAD \\", file = bashFile)
  print("  --family-expr 'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < 0.001' \\", file = bashFile)
  print("  --trio 'comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_popmax_af < 0.005' \\", file = bashFile)
  print("  | slivar_static compound-hets \\", file = bashFile)
  print("  --vcf /dev/stdin \\", file = bashFile)
  print("  --skip NONE \\", file = bashFile)
  print("  -s comphet_side \\", file = bashFile)
  print("  -s denovo \\", file = bashFile)
  print("  -p $PED \\", file = bashFile)
  print("  -o $SLIVAR2VCF \\", file = bashFile)
  print("  >> $STDOUT 2>> $STDERR", file = bashFile)
  print("  bcftools index -t $SLIVAR2VCF", file = bashFile)
  print(file = bashFile)

  # Delete the annotated vcf
  print("# Remove the annotated VCF", file = bashFile)
  print("  rm -f $ANNOTATEDVCF", file = bashFile)
  print("  rm -f \"$ANNOTATEDVCF\".tbi", file = bashFile)
  print(file = bashFile)

  # Combine the slivar vcf files
  versionLine = "##calypsoVersion=v" + version + "r" + resourceInfo['version']
  vcfLine     = "##calypsoSourceVcf=" + vcfFile
  print("# Combine the slivar vcf files", file = bashFile)
  print("  bcftools concat -a -d none $SLIVAR1VCF $SLIVAR2VCF \\", file = bashFile)
  print("  | bcftools annotate -H '" + versionLine + "' -H '" + vcfLine + "' -O z -o $FINALVCF \\", sep = "", file = bashFile)
  print("  >> $STDOUT 2>> $STDERR", file = bashFile)
  print("  bcftools index -t $FINALVCF", file = bashFile)
  print(file = bashFile)

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

  # Loop over all resources and see if any of them use VEP
  vepCommands = []
  vepFields   = []
  useVep      = False
  for resource in resourceInfo["resources"]:
    if resourceInfo["resources"][resource]["isVepAnnotation"] and resourceInfo["resources"][resource]["apply_to_filtered"]:
      useVep = True
      for command in resourceInfo["resources"][resource]["vep_commands"]: vepCommands.append(command)
      for field in resourceInfo["resources"][resource]["vep_fields"]: vepFields.append(field)

  # If so, include the VEP commands
  if useVep:
    print("  | vep \\", file = bashFile)
    print("    --vcf \\", file = bashFile)
    print("    --force \\", file = bashFile)
    print("    --check_existing \\", file = bashFile)
    print("    --coding_only \\", file = bashFile)
    print("    --quiet \\", file = bashFile)
    print("    --fork 40 \\", file = bashFile)
    print("    --format vcf \\", file = bashFile)
    print("    --force_overwrite \\", file = bashFile)
    print("    --cache \\", file = bashFile)
    print("    --no_stats \\", file = bashFile)
    print("    --fasta $REF \\", file = bashFile)
    print("    --dir_cache ", resourceInfo["resources"]["vep"]["cache"], " \\", sep = "", file = bashFile)
    print("    --dir_plugins ", resourceInfo["resources"]["vep"]["plugins"], " \\", sep = "", file = bashFile)
    print("    --assembly GRCh", str(args.reference), " \\", sep = "", file = bashFile)
    for command in vepCommands: print("    ", command, " \\", sep = "", file = bashFile)
    if len(vepFields) > 0: print("    --fields \"", ",".join(vepFields), "\" \\", sep = "", file = bashFile)
    print("    --output_file STDOUT \\", sep = "", file = bashFile)
    print("    2>> $STDERR \\", file = bashFile)
    print("  | bcftools view -O z -o $FILTEREDVCF - \\", file = bashFile)

  # If VEP is not used, finish the Slivar command
  else: print("  -o $FILTEREDVCF \\", file = bashFile)
  print("  >> $STDOUT 2>> $STDERR", file = bashFile)
  print(file = bashFile)

  # Extract vcf files of the de novo, x de novo, and recessive variants to add as variant sets. This can only be performed
  # for trios with a known proband
  if (args.family_type == "trio" or args.family_type == "quad") and proband: modeOfInheritance(vcfBase, bashFile)

  # Generate the tsv files to pass annotations to Mosaic
  print("# Generate the tsv files to pass annotations to Mosaic", file = bashFile)
  print("  echo \"Generating annotations tsv files for Mosaic\"", file = bashFile)
  print(file = bashFile)
  generateTsv(args, bashFile)

  # Close the file
  bashFile.close()

  # Make the annotation script executable
  makeExecutable = os.popen("chmod +x " + bashFilename).read()

  # Return the output vcf
  return finalVcf, filteredVcf

# Extract vcf files of the de novo, x de novo, and recessive variants to add as variant sets. This can only be performed
# for trios with a known proband
def modeOfInheritance(vcfBase, bashFile):
  global moiFiles

  moiFiles["denovoVcf"]    = {"file": str(vcfBase) + "_calypso_filtered_denovo.vcf.gz", "description": "Slivar de novo variants"}
  moiFiles["xdenovoVcf"]   = {"file": str(vcfBase) + "_calypso_filtered_xdenovo.vcf.gz", "description": "Slivar X de novo variants"}
  moiFiles["recessiveVcf"] = {"file": str(vcfBase) + "_calypso_filtered_recessive.vcf.gz", "description": "Slivar recessive variants"}
  #moiFiles["comphetVcf"]   = {"file": str(vcfBase) + "_calypso_filtered_comphet.vcf.gz", "description": "Slivar compound heterozygous variants"}
  print("# Generate vcf files containing variants of different modes of inheritance", file = bashFile)

  # Generate a vcf of de novo variants
  print("  # De novo variants", file = bashFile)
  print("  echo \"Generating vcf of de novo variants...\"", file = bashFile)
  print("  DENOVO_VCF=", moiFiles["denovoVcf"]["file"], sep = "", file = bashFile)
  print("  bcftools view -i 'INFO/denovo=\"", proband, "\"' -O z -o $DENOVO_VCF $FILTEREDVCF", sep = "", file = bashFile)
  print(file = bashFile)

  # Generate a vcf of X de novo variants
  print("  # X de novo variants", file = bashFile)
  print("  echo \"Generating vcf of x de novo variants...\"", file = bashFile)
  print("  X_DENOVO_VCF=", moiFiles["xdenovoVcf"]["file"], sep = "", file = bashFile)
  print("  bcftools view -i 'INFO/x_denovo=\"", proband, "\"' -O z -o $X_DENOVO_VCF $FILTEREDVCF", sep = "", file = bashFile)
  print(file = bashFile)

  # Generate a vcf of de novo variants
  print("  # De novo variants", file = bashFile)
  print("  echo \"Generating vcf of recessive variants...\"", file = bashFile)
  print("  RECESSIVE_VCF=", moiFiles["recessiveVcf"]["file"], sep = "", file = bashFile)
  print("  bcftools view -i 'INFO/recessive=\"", proband, "\"' -O z -o $RECESSIVE_VCF $FILTEREDVCF", sep = "", file = bashFile)
  print(file = bashFile)

################
################
################ NEED TO SORT COMP HETS
################
################
#  # Generate a vcf of comp het variants
#  print("  # Compound heterozygous variants", file = bashFile)
#  print("  echo \"Generating vcf of compound heterozygous variants...\"", file = bashFile)
#  print("  COMPHET_VCF=", moiFiles["comphetVcf"]["file"], sep = "", file = bashFile)
#  print("  bcftools view -i 'INFO/slivar_comphet=\"", proband, "\"' -O z -o $COMPHET_VCF $FILTEREDVCF", sep = "", file = bashFile)
#  print(file = bashFile)

# Generate the tsv files to pass annotations to Mosaic
def generateTsv(args, bashFile):
  global resourceInfo
  global mosaicInfo
  privateResources = []

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
      if annotation_type == "private": privateResources.append(resource)
      else: buildPublicTsv(args, resource, bashFile)

  # Build the tsv for private annotations
  if len(privateResources) > 0:

    # All of the private annotations need to be created in the project
    temp = 1
    annotationUids = createAnnotations(args, privateResources)
    buildPrivateAnnotationTsv(args, privateResources, annotationUids, bashFile)

# Build the tsv file
def buildPublicTsv(args, resource, bashFile):
  global mosaicInfo
  global tsvFiles

  # First create the header line
  header = "CHROM\\tSTART\\tEND\\tREF\\tALT"
  for annotation in sorted(mosaicInfo["resources"][resource]["annotations"]):
    header += "\\t" + mosaicInfo["resources"][resource]["annotations"][annotation]["uid"]

  # Print the header
  print("  # Public annotation: ", resource, sep = "", file = bashFile)
  outputFile         = resource + "_calypso_annotations_mosaic.tsv"
  tsvFiles[resource] = outputFile
  tempFile           = outputFile + ".tmp"
  print("  echo -e \"", header, "\" > ", outputFile, sep = "", file = bashFile)

  # Include all annotations for all resources
  annotations = "%CHROM\\t%POS\\t%END\\t%REF\\t%ALT"

  delimiter     = mosaicInfo["resources"][resource]["delimiter"]
  isPostprocess = mosaicInfo["resources"][resource]["postprocess"]

  # Loop over the annotations and add to the command line
  for annotation in sorted(mosaicInfo["resources"][resource]["annotations"]):
    if delimiter: annotations += "\\t%INFO/" + mosaicInfo["resources"][resource]["info"]
    else: annotations += "\\t%INFO/" + annotation

  # Build the query command
  print("  bcftools query -f '", annotations, "\\n' \\", sep = "", file = bashFile)
  print("  $FILTEREDVCF >> ", outputFile, sep = "", file = bashFile)
  print(file = bashFile)

  # If the resources need post processing (e.g. with additional scripts), do not modify the tsv
  if isPostprocess:
    print("  # Postprocess the annotation tsv", file = bashFile)
    postCommand = mosaicInfo["resources"][resource]["post_command"]
    command = (postCommand["precommand"] + " ") if postCommand["precommand"] else ""
    command += os.path.dirname( __file__) + "/" + postCommand["tool"]
    for argument in postCommand["args"]:
      argValue = postCommand["args"][argument]

      # Process the different arguments for standard values
      if argValue == "TSV": argValue = outputFile
      elif argValue == "data_directory": argValue = args.data_directory
      elif argValue == "reference": argValue = args.reference
      command += " " + argument + " " + argValue
    print("  ", command, sep = "", file = bashFile)
    print(file = bashFile)
  else:

    # If the delimiter has been set, the tsv contains the full annotation (e.g. A|B|C), for each
    # annotation. This needs to be modified to break the annotation up.
    print("  # Update the tsv file to Mosaic specifications", file = bashFile)
#    print("  awk '{FS=\"\\t\"; OFS=\"\\t\"} {if ($1 != \"CHROM\") {if ($1 ~ \"chr\") {$1=substr($1, 4)}; $3 = $3+1; " + awkCommand + "for(i=6; i<=NF; i++) {$i=($i==\".\" ? \"\":$i)}}; print $0}' \\", file = bashFile)
    print("  awk '{FS=\"\\t\"; OFS=\"\\t\"} {if ($1 != \"CHROM\") {if ($1 ~ \"chr\") {$1=substr($1, 4)}; $3 = $3+1; for(i=6; i<=NF; i++) {$i=($i==\".\" ? \"\":$i)}}; print $0}' \\", file = bashFile)
    print(" ", outputFile, ">", tempFile, file = bashFile)
    print("  mv", tempFile, outputFile, file = bashFile)
    print(file = bashFile)

# Creat new annotations
def createAnnotations(args, privateResources):
  global mosaicConfig
  global mosaicInfo
  global samples
  global sampleOrder
  global createdAnnotations
  annotationUids = {}

  # Get all the annotations in the project
  command = api_va.getVariantAnnotations(mosaicConfig, args.project_id, 100, 1)
  data    = json.loads(os.popen(command).read())

  # Store the annotations in the project by name
  projectAnnotations = {}
  for annotation in data: projectAnnotations[str(annotation["name"])] = annotation

  # Loop over the private resources
  for resource in sorted(privateResources):
    annotationUids[resource] = {}

    # Loop over all the annotations for this resource
    for annotation in mosaicInfo["resources"][resource]["annotations"]:

      # If this is a genotype annotation, create an attribute for each sample
      isGenotype = mosaicInfo["resources"][resource]["annotations"][annotation]["isGenotype"]
      annotationUids[resource][annotation] = {"isGenotype": isGenotype, "uids": []}
      if isGenotype:

        # Loop over all samples in the project and get the annotations for each sample
        for sample in sampleOrder:
          annotationName = str(annotation) + " " + str(samples[sample]["relationship"])

          # Check if an annotation with this name already exists
          annotationUid = projectAnnotations[annotationName]["uid"] if annotationName in projectAnnotations else False

          # If the annotation does not exist, create it
          if not annotationUid:
            valueType = mosaicInfo["resources"][resource]["annotations"][annotation]["type"]
            data      = json.loads(os.popen(api_va.postCreateVariantAnnotations(mosaicConfig, annotationName, valueType, "private", args.project_id)).read())

            # Get the uid of the newly created annotation
            annotationUid = data["uid"]

          # Store the annotation uids
          annotationUids[resource][annotation]["uids"].append(annotationUid)
          createdAnnotations[annotationName] = annotationUid

      # Other types of private annotations are not yet handled
      else: fail("PRIVATE ANNOTAION TYPE NOT HANDLED")

  return annotationUids

# Build the tsv containing the private annotations
def buildPrivateAnnotationTsv(args, resources, annotationUids, bashFile):
  global mosaicConfig
  global tsvFiles
  global sampleOrder

  # First create the header line
  header = "CHROM\\tSTART\\tEND\\tREF\\tALT"
  for resource in sorted(resources):
    for annotation in sorted(annotationUids[resource]):
      for uid in annotationUids[resource][annotation]["uids"]: header += "\\t" + uid

  # Print the header
  print("  # Private annotation: ", resource, sep = "", file = bashFile)
  outputFile         = resource + "_calypso_annotations_mosaic.tsv"
  tsvFiles[resource] = outputFile
  tempFile           = outputFile + ".tmp"
  print("  echo -e \"", header, "\" > ", outputFile, sep = "", file = bashFile)

  # Include all annotations for all resources
  annotations = "%CHROM\\t%POS\\t%END\\t%REF\\t%ALT"

  delimiter     = mosaicInfo["resources"][resource]["delimiter"]
  isPostprocess = mosaicInfo["resources"][resource]["postprocess"]

  # Loop over the annotations and add to the command line
  for annotation in sorted(mosaicInfo["resources"][resource]["annotations"]):

    # Determine if this is a genotype annotation, loop over the samples in the correct order, and add the annotation
    if mosaicInfo["resources"][resource]["annotations"][annotation]["isGenotype"]: annotations += "[\\t%GQ]"
    elif delimiter: annotations += "\\t%INFO/" + mosaicInfo["resources"][resource]["info"]
    else: annotations += "\\t%INFO/" + annotation

  # Build the query command
  print("  bcftools query -f '", annotations, "\\n' \\", sep = "", file = bashFile)
  print("  $FILTEREDVCF >> ", outputFile, sep = "", file = bashFile)
  print(file = bashFile)

  # If the resources need post processing (e.g. with additional scripts), do not modify the tsv
  if isPostprocess:
    print("  # Postprocess the annotation tsv", file = bashFile)
    postCommand = mosaicInfo["resources"][resource]["post_command"]
    command = (postCommand["precommand"] + " ") if postCommand["precommand"] else ""
    command += os.path.dirname( __file__) + "/" + postCommand["tool"]
    for argument in postCommand["args"]:
      argValue = postCommand["args"][argument]

      # Process the different arguments for standard values
      if argValue == "TSV": argValue = outputFile
      elif argValue == "data_directory": argValue = args.data_directory
      elif argValue == "reference": argValue = args.reference
      command += " " + argument + " " + argValue
    print("  ", command, sep = "", file = bashFile)
    print(file = bashFile)
  else:

    # If the delimiter has been set, the tsv contains the full annotation (e.g. A|B|C), for each
    # annotation. This needs to be modified to break the annotation up.
    print("  # Update the tsv file to Mosaic specifications", file = bashFile)
#    print("  awk '{FS=\"\\t\"; OFS=\"\\t\"} {if ($1 != \"CHROM\") {if ($1 ~ \"chr\") {$1=substr($1, 4)}; $3 = $3+1; " + awkCommand + "for(i=6; i<=NF; i++) {$i=($i==\".\" ? \"\":$i)}}; print $0}' \\", file = bashFile)
    print("  awk '{FS=\"\\t\"; OFS=\"\\t\"} {if ($1 != \"CHROM\") {if ($1 ~ \"chr\") {$1=substr($1, 4)}; $3 = $3+1; for(i=6; i<=NF; i++) {$i=($i==\".\" ? \"\":$i)}}; print $0}' \\", file = bashFile)
    print(" ", outputFile, ">", tempFile, file = bashFile)
    print("  mv", tempFile, outputFile, file = bashFile)
    print(file = bashFile)

# Output a script to upload variants to Mosaic
def uploadVariants(args, filteredVcf):
  global mosaicInfo
  global workingDir
  global proband
  global date

  # Open a script file
  uploadFileName  = workingDir + "calypso_upload_variants_to_mosaic.sh"
  try: uploadFile = open(uploadFileName, "w")
  except: fail("Could not open " + str(uploadFileName) + " to write to")

  # Define the name of the variant set that will be create
  description = "Calypso_v" + str(mosaicInfo["version"]) + "_variants_" + date

  # Write the command to file
  print("# Upload variants to Mosaic", file = uploadFile)
  print("python ", os.path.dirname( __file__ ), "/mosaic_commands/upload_variants_to_mosaic.py \\", sep = "", file = uploadFile)
  print("  -i ", str(filteredVcf), " \\", sep = "", file = uploadFile)
  print("  -a ", os.path.dirname( __file__), "/mosaic_commands \\", sep = "", file = uploadFile)
  print("  -d \"", description, "\" \\", sep = "", file = uploadFile)
  if args.config: print("  -c", str(args.config), "\\", file = uploadFile)
  else: print("  -c \"Insert config file here\" \\", sep = "", file = uploadFile)
  if args.project_id: print("  -p", str(args.project_id), file = uploadFile)
  else: print("  -p \"Insert Mosaic project id here\"", file = uploadFile)
  print(file = uploadFile)

  # If vcf files of different modes of inheritance were generated, upload these too
  if (args.family_type == "trio" or args.family_type == "quad") and proband:
    print("# Upload the mode of inheritance files as variant sets", file = uploadFile)
    for moiFile in moiFiles:
      print("python ", os.path.dirname( __file__ ), "/mosaic_commands/upload_variants_to_mosaic.py \\", sep = "", file = uploadFile)
      print("  -i ", str(moiFiles[moiFile]["file"]), " \\", sep = "", file = uploadFile)
      print("  -a ", os.path.dirname( __file__), "/mosaic_commands \\", sep = "", file = uploadFile)
      print("  -d \"", moiFiles[moiFile]["description"], "\" \\", sep = "", file = uploadFile)
      if args.config: print("  -c", str(args.config), "\\", file = uploadFile)
      else: print("  -c \"Insert config file here\" \\", sep = "", file = uploadFile)
      if args.project_id: print("  -p", str(args.project_id), file = uploadFile)
      else: print("  -p \"Insert Mosaic project id here\"", file = uploadFile)
      print(file = uploadFile)

  # Close the file
  uploadFile.close()

  # Make the annotation script executable
  makeExecutable = os.popen("chmod +x " + str(uploadFileName)).read()

# Output scripts to upload annotations to Mosaic
def uploadAnnotations(args):
  global tsvFiles
  global mosaicInfo
  global workingDir
  uploadFileName  = workingDir + "calypso_upload_annotations_to_mosaic.sh"
  isScript        = False

  # Loop over all resources
  for resource in tsvFiles:
    annotationType = mosaicInfo["resources"][resource]["annotation_type"]
    isComplete     = False

    # Public annotations need to be imported into the Mosaic project
    uids = {}
    if annotationType == "public": isComplete = getPublicAnnotation(args, resource)
    elif annotationType == "private": isComplete = True

    # Only create the script file if the annotation(s) exists in Mosaic
    if isComplete: 

      # Create a single script file to upload all variants
      if not isScript: 
        try: uploadFile = open(uploadFileName, "w")
        except: fail("Could not open " + str(uploadFileName) + " to write to")
        isScript = True

      # Write the command to file
      print("# Upload ", resource, " annotations to Mosaic", file = uploadFile)
      print("python ", os.path.dirname( __file__ ), "/mosaic_commands/upload_annotations_to_mosaic.py \\", sep = "", file = uploadFile)
      print("  -i ", str(tsvFiles[resource]), " \\", sep = "", file = uploadFile)
      print("  -a ", os.path.dirname( __file__), "/mosaic_commands \\", sep = "", file = uploadFile)
      if args.config: print("  -c", str(args.config), "\\", file = uploadFile)
      else: print("  -c \"Insert config file here\" \\", sep = "", file = uploadFile)
  
      # Public annotations need to be uploaded to the project that contains them, private
      # annotations are uploaded to the defined project.
      if annotationType == "private":
        if args.project_id: print("  -p", str(args.project_id), file = uploadFile)
        else: print("  -p \"Insert Mosaic project id here\"", file = uploadFile)
      else: print("  -p ", mosaicInfo["resources"][resource]["project_id"], sep = "", file = uploadFile)
      print(file = uploadFile)
  
  # If a script files was created, close it and make it executable
  if isScript:

    # Close the file
    uploadFile.close()
  
    # Make the annotation script executable
    makeExecutable = os.popen("chmod +x " + str(uploadFileName)).read()

# Check the public annotations exist, and get their ids
def getPublicAnnotation(args, resource):
  global mosaicInfo
  global warnings
  page = 1
  uids = {}

  # Find the annotation id for the public annotation. Get the first page of annotations and check
  # if there are additional pages required
  for annotation in mosaicInfo["resources"][resource]["annotations"]:
    uids[mosaicInfo["resources"][resource]["annotations"][annotation]["uid"]] = False

  # Get the first page of public attributes, and how many pages there are
  noPages          = getAnnotationsPages(args)
  isComplete, uids = getAnnotations(args, page, uids)

  # If some ids are not found, go to the next page. If there are no more pages, throw a warning
  while page <= noPages:
    isComplete, uids = getAnnotations(args, page, uids)
    page += 1
    if isComplete: break

  # If the annotations weren't found, throw a warning
  if not isComplete:
    warningTitle = "Annotation not found"
    description  = "Public annotations need to be imported into the Mosaic project, but these annotations weren't found "
    description += "and so cannot be imported"
    if warningTitle not in warnings: warnings[warningTitle] = {"description": description, "messages": []}
    for annotation in uids:
      if not uids[annotation]: warnings[warningTitle]["messages"].append("  " + annotation)

  # If the annotations can be imported, import them
  else:
    for annotation in uids: importAnnotation(args, annotation, uids[annotation])

  # Return the isComplete variable
  return isComplete

# Determine how many pages of public attributes there are
def getAnnotationsPages(args):
  global mosaicConfig

  # Get a page of public annotations from Mosaic (100 annotations per page)
  command = api_va.getVariantAnnotationsImport(mosaicConfig, args.project_id, 100, 1)
  data    = json.loads(os.popen(command).read())
  noPages = math.ceil(int(data["count"]) / int(100))

  # Return the number of pages of annotations
  return noPages

# Get a page of public attributes
def getAnnotations(args, page, uids):
  global mosaicConfig

  # Get a page of public annotations from Mosaic (100 annotations per page)
  command = api_va.getVariantAnnotationsImport(mosaicConfig, args.project_id, 100, page)
  data    = json.loads(os.popen(command).read())

  for annotation in data["data"]:
    if annotation["uid"] in uids.keys(): uids[annotation["uid"]] = annotation["id"]

  # Check if all required attributes have been found and the id stored
  isComplete = True
  for annotation in uids:
    if not uids[annotation]: isComplete = False

  # Return the Mosaic uids and whether they were all found
  return isComplete, uids

# Import a public annotation into the Mosaic project
def importAnnotation(args, annotation, uid):
  global mosaicConfig

  # Create the command to import the annotation
  command = api_va.postImportVariantAnnotations(mosaicConfig, uid, args.project_id)
  data    = json.loads(os.popen(command).read())

# Output summary file
def calypsoSummary(args, finalVcf):
  global resourceInfo
  global workingDir
  global version
  global date

  # Open an output summary file
  summaryFileName  = workingDir + "calypso_" + date + ".txt"
  try: summaryFile = open(summaryFileName, "w")
  except: fail("Failed to open summary file")

  # Write relevant information to file
  print("### Output from Calypso pipeline ###", file = summaryFile)
  print(file = summaryFile)
  print("Calypso pipeline version: ", version, sep = "", file = summaryFile)
  print("Calypso resource version: ", resourceInfo['version'], sep = "", file = summaryFile)
  print("Reference:                ", args.reference, sep = "", file = summaryFile)
  print("Created on:               ", date, sep = "", file = summaryFile)
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

# Set project attributes to indicate when and which version of Calypso has been run
def updateCalypsoAttributes(args):
  global mosaicConfig
  global version
  global date

  # Define the public attributes we want to import or update
  calypso = {}
  calypso["Calypso version"]  = {"uid": False, "id": False, "value": version, "inProject": False}
  calypso["Calypso date run"] = {"uid": False, "id": False, "value": date, "inProject": False}
  calypso["Calypso history"]  = {"uid": False, "id": False, "value": str(version) + ":" + date, "inProject": False}

  # Get all project attributes to check if the Calypso attributes have already been imported and set
  command = api_pa.getPublicProjectAttributes(mosaicConfig, 100, 1)
  data    = json.loads(os.popen(command).read())
  noPages = math.ceil(int(data["count"]) / int(100))

  # There may be multiple pages of attributes, so loop over all pages looking for the attributes we need
  page = 1
  while page <= noPages:
    command = api_pa.getPublicProjectAttributes(mosaicConfig, 100, 1)
    data    = json.loads(os.popen(command).read())
    page   += 1

    # Loop over the attributes and look for the Calypso attributes
    for attribute in data["data"]:
      if attribute["name"] in calypso: 
        calypso[attribute["name"]]["uid"]   = attribute["uid"]
        calypso[attribute["name"]]["id"]    = attribute["id"]

  # Loop over the attributes in the current project and check for the Calypso attributes
  command = api_pa.getProjectAttributes(mosaicConfig, args.project_id)
  data    = json.loads(os.popen(command).read())
  for attribute in data:

    # If an attribute with the right name exists, check the Mosaic id matches the public attribute, and
    # update the value
    if attribute["name"] in calypso:
      if attribute["id"] == calypso[attribute["name"]]["id"]:

        # If this is the history, the current version and date should be appended. If the last execution of the
        # Calypso pipeline was today and with the current version, however, do not update.
        if attribute["name"] == "Calypso history": calypso[attribute["name"]]["value"] = updateHistory(args, attribute["values"])
        command = api_pa.putProjectAttribute(mosaicConfig, calypso[attribute["name"]]["value"], args.project_id, attribute["id"])
        data    = json.loads(os.popen(command).read())
        calypso[attribute["name"]]["inProject"] = True

  # Import all attributes that are not in the project
  for attribute in calypso:
    if not calypso[attribute]["inProject"]:
      command = api_pa.postImportProjectAttribute(mosaicConfig, calypso[attribute]["id"], calypso[attribute]["value"], args.project_id)
      data    = json.loads(os.popen(command).read())

# Update the Calypso history value
def updateHistory(args, valueRecords):
  global version
  global date

  # Get the history string for this project
  for record in valueRecords:
    value = record["value"] if str(record["project_id"]) == str(args.project_id) else False

  # If there was no history value, return todays date and version
  if not value: return str(version) + ":" + date

  # Get the most recent update to the history
  mostRecent  = value.split(",")[-1] if "," in value else value
  lastVersion = mostRecent.split(":")[0]
  lastDate    = mostRecent.split(":")[1]

  # If the date and version at execution are the same as at the last execution, do not update
  if str(lastVersion) == str(version) and str(lastDate) == date: return value
  else: return str(value) + "," + str(version) + ":" + date

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
version = "0.2.0"
date    = str(date.today())

# The working directory where all created files are kept
workingDir = os.getcwd() + "/calypso_v" + version + "r"

# Store warnings to be output at the end
warnings = {}

# Store info on allowed values
isFamily          = True
allowedFamily     = ['singleton', 'duo', 'trio', 'quad', 'quintet']
allowedReferences = ['37', '38']

# Store the proband id
proband = False

# Store the mode of inheritance file names
moiFiles = {}

# Store information about the resource files
resourceVersion = False
resourceInfo    = {}

# Store information related to Mosaic
mosaicConfig = {}
mosaicInfo   = {}

# Store information on the project samples
samples     = {}
sampleOrder = []

# Store the tsv files that upload annotations to Mosaic
tsvFiles = {}

# Store the path to the Calypso directory
rootPath = os.path.dirname( __file__)

# Store the allowed genotype options for saved filters
genotypeOptions = []
genotypeOptions.append("ref_samples")
genotypeOptions.append("alt_samples")
genotypeOptions.append("het_samples")
genotypeOptions.append("hom_samples")

# Store annotations created for this project
createdAnnotations = {}

# Created files
tomlFilename = ""

if __name__ == "__main__":
  main()
