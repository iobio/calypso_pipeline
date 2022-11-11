#!/usr/bin/python

from __future__ import print_function
from datetime import date
from os.path import exists
from sys import path

import os
import argparse
import json
import math
import glob
import importlib

# Add the path of the common functions and import them
path.append(os.path.dirname(__file__) + "/components")
import calypso_bash_script as bashScript
import calypso_mosaic as mos
import calypso_mosaic_resources as mosr
import calypso_resources as res
import calypso_samples as sam
import calypso_toml as toml
import calypso_upload as upload
import calypso_variant_filters as vfilt

def main():
  global mosaicConfig
  global workingDir

  # Parse the command line
  args = parseCommandLine()

  # Check the supplied parameters are as expected, then expand the path to enable scripts and api commands
  # to be accessed for Calypso
  resourceInfo = checkArguments(args)
  buildPath(args)

  # Read the resources json file to identify all the resources that will be used in the annotation pipeline
  resourceInfo = res.readResources(args, resourceInfo)
  rootPath = os.path.dirname( __file__)
  setWorkingDir(resourceInfo['version'])

  # Import additional tools
  import mosaic_config
  import api_project_attributes as api_pa
  import api_project_settings as api_ps
  import api_samples as api_s
  import api_variants as api_v
  import api_variant_annotations as api_va
  import api_variant_filters as api_vf

  # Parse the Mosaic config file to get the token and url for the api calls
  mosaicRequired = {"token": True, "url": True, "attributesProjectId": True}
  mosaicConfig   = mosaic_config.parseConfig(args, mosaicRequired)

  # Determine the id of the proband, and the order of the samples in the vcf header
  proband, samples = sam.getProband(mosaicConfig, args.ped, args.family_type, args.project_id, api_s)
  sam.getSampleOrder(samples, args.input_vcf)

  # Get all the project attributes in the target Mosaic project
  projectAttributes = mos.getProjectAttributes(mosaicConfig, args.project_id, api_pa)
  publicAttributes  = mos.getPublicProjectAttributes(mosaicConfig, api_pa)

  # Read the Mosaic json and validate its contents
  mosaicInfo         = mosr.readMosaicJson(args)
  projectAnnotations = mos.getProjectAnnotations(mosaicConfig, args.project_id, api_va)
  publicAnnotations  = mos.getPublicAnnotations(mosaicConfig, args.project_id, api_va)
  privateAnnotations = mos.createPrivateAnnotations(mosaicConfig, mosaicInfo['resources'], projectAnnotations, samples, args.project_id, api_va)
  mosr.checkPublicAnnotations(mosaicConfig, mosaicInfo['resources'], publicAnnotations, args.project_id, api_va)

  # Build the toml file for vcfanno. This defines all of the annotations that vcfanno needs to use
  # with the path to the required files
  tomlFilename = toml.buildToml(workingDir, resourceInfo)

  # Generate the bash script to run the annotation pipeline
  bashFilename, bashFile = bashScript.openBashScript(workingDir)
  filteredVcf = bashScript.bashResources(resourceInfo, workingDir, bashFile, args.input_vcf, args.ped, tomlFilename)
  bashScript.samplesFile(bashFile)
  bashScript.cleanedVcf(bashFile)
  bashScript.annotateVcf(bashFile, resourceInfo)
  bashScript.filterVariants(bashFile, proband, resourceInfo)
  bashScript.rareDiseaseVariants(bashFile)
  bashScript.deleteFiles(bashFile)
  bashScript.finishScript(bashFile, bashFilename)

  # Prepare the target project. This includes:
  # 1. Remove any unnecessary annotations from the project
  # 2. Import public annotations into the project
  # 3. Set the default annotations for the project
  # 4. Set up required variant filters
  # 5. Add and set project attributes
  if 'remove' in mosaicInfo: mos.removeAnnotations(mosaicConfig, mosaicInfo['remove'], projectAnnotations, args.project_id, api_va)
  mos.importAnnotations(mosaicConfig, mosaicInfo['resources'], projectAnnotations, publicAnnotations, args.project_id, api_va)
  mos.defaultAnnotations(mosaicConfig, mosaicInfo['defaultAnnotations'], publicAnnotations, privateAnnotations, args.project_id, api_ps)
  vfilt.variantFilters(mosaicConfig, rootPath, samples, privateAnnotations, args.project_id, api_vf)
  mos.updateCalypsoAttributes(mosaicConfig, projectAttributes, publicAttributes, version, args.project_id, api_pa)

  # Generate scripts to upload filtered variants and annotations to Mosaic
  upload.uploadVariants(mosaicInfo, mosaicConfig, workingDir, args.config, filteredVcf, args.project_id, api_v)
  #upload.uploadAnnotations(args)

#  # If a previous Calypso vcf file was set, get the version assigned to the vcf
#  if args.previous_vcf: getPreviousVcfInfo(args)

  # Output a summary file listing the actions undertaken by Calypso with all version histories
#  calypsoSummary(args, finalVcf)

# Input options
def parseCommandLine():
  global version
  global allowedReferences
  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  parser.add_argument('--reference', '-r', required = True, metavar = 'string', help = 'The reference genome to use. Allowed values: ' + ', '.join(allowedReferences))
  parser.add_argument('--family_type', '-f', required = True, metavar = 'string', help = 'The familty structure. Allowed values: ' + ', '. join(allowedFamily))
  parser.add_argument('--data_directory', '-d', required = True, metavar = 'string', help = 'The path to the directory where the resources live')
  parser.add_argument('--utils_directory', '-l', required = True, metavar = 'string', help = 'The path to the public-utils directory')
  parser.add_argument('--input_vcf', '-i', required = True, metavar = 'string', help = 'The input vcf file to annotate')
  parser.add_argument('--ped', '-e', required = True, metavar = 'string', help = 'The pedigree file for the family. Not required for singletons')
  parser.add_argument('--project_id', '-p', required = True, metavar = 'string', help = 'The project id that variants will be uploaded to')
  parser.add_argument('--config', '-c', required = True, metavar = "string", help = "The config file for Mosaic")

  # Optional reanalysis arguments
#  parser.add_argument('--previous_vcf', '-s', required = False, metavar = "string", help = "A previously annotated Calypso vcf file")

  # Optional pipeline arguments
  parser.add_argument('--resource_json', '-j', required = False, metavar = 'string', help = 'The json file describing the annotation resources')
#  parser.add_argument('--no_comp_het', '-n', required = False, action = "store_true", help = "If set, comp het determination will be skipped")
#
  # Optional argument to handle HPO terms
#  parser.add_argument('--hpo', '-o', required = False, metavar = "string", help = "A comma separate list of hpo ids for the proband")

  # Optional mosaic arguments
  parser.add_argument('--token', '-t', required = False, metavar = "string", help = "The Mosaic authorization token")
  parser.add_argument('--url', '-u', required = False, metavar = "string", help = "The base url for Mosaic")
  parser.add_argument('--attributes_project', '-a', required = False, metavar = "integer", help = "The Mosaic project id that contains public attributes")
  parser.add_argument('--mosaic_json', '-m', required = False, metavar = 'string', help = 'The json file describing the Mosaic parameters')

  # Version
  parser.add_argument('--version', '-v', action="version", version='Calypso annotation pipeline version: ' + str(version))

  return parser.parse_args()

# Check the supplied arguments
def checkArguments(args):
  global allowedReferences
  global allowedFamily
  resourceInfo = {}

  # Check the reference is allowed
  if args.reference not in allowedReferences:
    message = 'The allowed references (--reference, -r) are:'
    for ref in allowedReferences: message += '\n  ' + str(ref)
    fail(message)

  # Ensure the data path ends with a "/", then add the reference directory
  if args.data_directory[-1] != "/": args.data_directory += "/"
  resourceInfo['path'] = args.data_directory + "GRCh" + str(args.reference) + "/"

  # Check that the public-utils directory exists and add to the path so scripts from here can be used
  if args.utils_directory[-1] != "/": args.utils_directory += "/"
  if not os.path.exists(args.utils_directory): fail('public-utils directory could not be found')

  # Define the name of the resource description json file. Use the file provided on the command
  # line, and resort to the default file if not included. The version of the resource file is
  # unknown, so look for any json with the correct prefix and if there is more than one such
  # file, throw an error.
  if args.resource_json: resourceFilename = args.resource_json
  else:
    resourceFilename = args.data_directory + 'resources_GRCh' + str(args.reference)
    resourceFiles    = glob.glob(resourceFilename + "*json")
    if len(resourceFiles) != 1: fail("There are zero, or more than one resource files for GRCh" + args.reference + " in " + args.data_directory)
    resourceFilename = resourceFiles[0]
  resourceInfo['json'] = resourceFilename

  # Check the family type is allowed
  if args.family_type not in allowedFamily: fail('The allowed family types (--family-type, -f) are: ' + ', '.join(allowedFamily))

  # If this is a singleton, no family based analyis should be performed
  resourceInfo['family_type'] = args.family_type
  resourceInfo['isFamily'] = False if args.family_type == 'singleton' else True

  # Check that the input files exist and have the correct extensions
  if not exists(args.input_vcf): fail('The vcf file could not be found')
  elif not args.input_vcf.endswith('.vcf.gz'): fail('The input vcf file (--vcf, -v) must be a compressed, indexed vcf and have the extension ".vcf.gz"')
  if resourceInfo['isFamily']:
    if not args.ped: fail('A ped file needs to specified (--ped, -e) for family type "' + str(args.family_type) + '"')
    if not exists(args.ped): fail('The ped file could not be found')
    elif not args.ped.endswith('.ped'): fail('The input ped file (--ped, -p) must have the extension ".ped"')

  # Return the info on resources
  return resourceInfo

# Build the path to allow importing additional modules
def buildPath(args):

  # Add public-utils and required directories to the path
  path.append(args.utils_directory)
  path.append(args.utils_directory + 'common_components')
  path.append(args.utils_directory + 'scripts')
  path.append(args.utils_directory + 'api_commands')

# Create a directory where all Calypso associated files will be stored
def setWorkingDir(version):
  global workingDir

  workingDir += version + "/"

  # If the directory doesn't exist, create it
  if not os.path.exists(workingDir): os.makedirs(workingDir)






















  
## If a previous Calypso vcf file was set, get the version assigned to the vcf
#def getPreviousVcfInfo(args):
#  global resourceInfo
#  global workingDir
#  global version
#  headerVcf = workingDir + "calypso_previous_header.vcf"
#  previousVersions = []
#  previousFiles    = []
#
#  # Check that the previous vcf file exists
#  if not exists(args.previous_vcf): fail("Previous Calypso vcf file, " + args.previous_vcf + ", does not exist.")
#
#  # Grab the vcf header and extract the previous versions
#  try: os.popen("bcftools view -h -O v -o " + headerVcf + " " + args.previous_vcf).read()
#  except: fail("Unable to extract header from vcf file, " + args.previous_vcf)
#  try: headerInfo = open(headerVcf, "r")
#  except: fail("Unable to read header information from file " + arg.previous_vcf)
#  for line in  headerInfo.readlines():
#    if line.rstrip().startswith("##calypsoVersion="): previousVersions.append(line.rstrip().split("=")[1])
#    if line.rstrip().startswith("##calypsoSourceVcf="): previousFiles.append(line.rstrip().split("=")[1])
#
#  # If there are multiple Calypso versions in the header, instruct the user to clean this up to only include the most recent
#  if len(previousVersions) != 1: fail("The previous Calypso vcf file, " + args.previous_vcf + ", either has no, or multiple \"##calypsoVersion=VERSION\" header lines. Ensure it only has one")
#  if len(previousFiles) != 1: fail("The previous Calypso vcf file, " + args.previous_vcf + ", either has no, or multiple \"##calypsoSourceVcf=VCFFILE\" header lines. Ensure it only has one")
#
#  # Get the previous Calypso and the resource version
#  versions = previousVersions[0].split("r")
#  resourceInfo['previous_version'] = versions[1]
#  previousVersion = versions[0].strip('v')
#  previousFile    = previousFiles[0].split('/')[-1]
#
#  # If both the Calypso and the resource versions are the same as for the previously generated vcf, and the source
#  # vcf files are the same, no changes have been made, so there is no need to rerun the pipeline
#  if str(resourceInfo['version']) == str(resourceInfo['previous_version']) and str(previousVersion) == str(version) and str(args.input_vcf) == str(previousFile):
#    failString  = "Running the same Calypso version (" + version + "), with the same resource version (" + resourceInfo['version'] + "), on the same\nsource vcf ("
#    failString += args.input_vcf + ") will result in the same output, so Calypso will not be run"
#    fail(failString)
#
#  # Delete the created header file
#  try: os.remove(headerVcf)
#  except: fail("Failed to delete created header vcf, " + headerVcf)
#
#
## Create variant filters in Mosaic
#def variantFilters(args):
#  global mosaicConfig
#  global samples
#  global rootPath
#  global genotypeOptions
#
#  # Get information about the samples
#  probandId = False
#  motherId  = False
#  fatherId  = False
#  jsonPath  = rootPath + str("/mosaic_filters/")
#
#  # Get the id of the proband and the parents
#  for sample in samples:
#    if samples[sample]["relationship"] == "Proband": probandId = samples[sample]["mosaicId"]
#    elif samples[sample]["relationship"] == "Mother": motherId = samples[sample]["mosaicId"]
#    elif samples[sample]["relationship"] == "Father": fatherId = samples[sample]["mosaicId"]
#
#  # Get all of the filters that exist in the project
#  existingFilters = {}
#  try:
#    for existingFilter in json.loads(os.popen(api_vf.getVariantFilters(mosaicConfig, args.project_id)).read()): existingFilters[existingFilter["name"]] = existingFilter["id"]
#  except: fail("Unable to get existing variant filters for project " + str(args.project_id))
#  
#  # Loop over all the filter json files
#  for filterJson in os.listdir(jsonPath):
#    if filterJson.endswith(".json"):
#
#      try: filterFile = open(jsonPath + filterJson, "r")
#      except: fail("File, " + jsonPath + filterJson + ", could not be opened")
#      try: data = json.load(filterFile)
#      except: fail("File, " + jsonPath + filterJson + ", is not a valid json file")
#      filterFile.close()
#
#      # Get the name of the filter to create
#      try: filterName = data["name"]
#      except: fail("Varianr filter file, " + str(filterJson) + ", defines a filter with no name. A name needs to be supplied")
#
#      # Check if the filter already exists in the project. If so, remove it
#      if filterName in existingFilters:
#        try: deleteFilter = os.popen(api_vf.deleteVariantFilter(mosaicConfig, args.project_id, existingFilters[filterName])).read()
#        except: fail("Unable to delete variant filter, " + str(filterName) + " from project " + str(args.project_id))
#
#      # Check what genotype filters need to be applied
#      try: genotypes = data["genotypes"]
#      except: fail("Mosaic variant filter json file, " + str(filterJson) + ", does not contain the \"genotypes\" field")
#      for geno in genotypes:
#        if geno not in genotypeOptions: fail("Mosaic variant filter json file, " + str(filterJson) + ", contains unknown genotype options: " + str(geno))
#        if not genotypes[geno]: continue
#
#        # Check which samples need to have the requested genotype and add to the command
#        sampleList = []
#        for sample in genotypes[geno]:
#          if sample == "proband": sampleList.append(probandId)
#          elif sample == "mother":
#            if motherId: sampleList.append(motherId)
#          elif sample == "father":
#            if fatherId: sampleList.append(fatherId)
#          elif sample == "sibling":
#            fail("Can't handle siblings in filters yet")
#          else: fail("Unknown sample, " + str(sample) + " in genotypes for Mosaic variant filter (" + str(filterJson) + ")")
#
#        # Add the genotype filter to the filters listed in the json
#        data["filters"][geno] = sampleList
#
#      # Check the filters provided in the json. The annotation filters need to extract the
#      # uids for the annotations
#      removeAnnotations = []
#      for index, annotationFilter in enumerate(data["filters"]["annotation_filters"]):
#
####################
####################
#################### ADD CHECKS ON THE FILTER JSON
####################
####################
#
#        # If the uid is provided with the filter, this annotation is complete. If not, the name
#        # field must be present, and this must point to a created annotation uid. Remove the name
#        # field and replace it with the uid.
#        if "uid" not in annotationFilter:
#          try: annotationName = annotationFilter["name"]
#          except: fail("Annotation " + str(filterJson) + " does not have a uid, but also no name")
#          annotationFilter.pop("name")
#
#          # If the annotation is listed as optional and the name doesn't point to a created annotation
#          # delete this annotation from the filter. This could be, for example, a filter on genotype
#          # quality for the mother, but the mother isn't present in the project
#          isOptional = annotationFilter.pop("optional") if "optional" in annotationFilter else False
#          try: annotationFilter["uid"] = createdAnnotations[annotationName]
#          except: 
#            if isOptional: removeAnnotations.append(index)
#            else: fail("Annotation, " + str(annotationName) + ", in " + str(filterJson) + " does not have a uid, and the name does not point to a created annotation")
#
#      # Remove any optional filters that did not have uids
#      for index in reversed(removeAnnotations): data["filters"]["annotation_filters"].pop(index)
#
#      # Post a new filter
#      execute = json.loads(os.popen(api_vf.postVariantFilter(mosaicConfig, data["name"], data["filters"], args.project_id)).read())
#
##################
##################
################## Do we want to extract all mois and upload as variant sets?
##################
##################
#  # Extract vcf files of the de novo, x de novo, and recessive variants to add as variant sets. This can only be performed
#  # for trios with a known proband
##  if (args.family_type == "trio" or args.family_type == "quad") and proband: modeOfInheritance(vcfBase, bashFile)
#
#  # Generate the tsv files to pass annotations to Mosaic
#  print("# Generate the tsv files to pass annotations to Mosaic", file = bashFile)
#  print("  echo \"Generating annotations tsv files for Mosaic\"", file = bashFile)
#  print(file = bashFile)
#  generateTsv(args, bashFile)
#
#  # If HPO ids were supplied, add a command to generate hpo annotations
#  if args.hpo: processHpo(args, bashFile)
#
#  # Closing notification
#  print("  echo \"Everything completed! Annotated VCF written to $FINALVCF\"", file = bashFile)
#  print(file = bashFile)
#
#  # Close the file
#  bashFile.close()
#
#  # Make the annotation script executable
#  makeExecutable = os.popen("chmod +x " + bashFilename).read()
#
#  # Return the output vcf
#  return finalVcf, filteredVcf
#
## Extract vcf files of the de novo, x de novo, and recessive variants to add as variant sets. This can only be performed
## for trios with a known proband
#def modeOfInheritance(vcfBase, bashFile):
#  global moiFiles
#
#  moiFiles["denovoVcf"]    = {"file": str(vcfBase) + "_calypso_filtered_denovo.vcf.gz", "description": "Slivar de novo variants"}
#  moiFiles["xdenovoVcf"]   = {"file": str(vcfBase) + "_calypso_filtered_xdenovo.vcf.gz", "description": "Slivar X de novo variants"}
#  moiFiles["recessiveVcf"] = {"file": str(vcfBase) + "_calypso_filtered_recessive.vcf.gz", "description": "Slivar recessive variants"}
#  moiFiles["comphetVcf"]   = {"file": str(vcfBase) + "_calypso_filtered_comphet.vcf.gz", "description": "Slivar compound heterozygous variants"}
#  print("# Generate vcf files containing variants of different modes of inheritance", file = bashFile)
#
#  # Generate a vcf of de novo variants
#  print("  # De novo variants", file = bashFile)
#  print("  echo \"Generating vcf of de novo variants...\"", file = bashFile)
#  print("  DENOVO_VCF=", moiFiles["denovoVcf"]["file"], sep = "", file = bashFile)
#  print("  bcftools view -i 'INFO/denovo=\"", proband, "\"' -O z -o $DENOVO_VCF $FILTEREDVCF", sep = "", file = bashFile)
#  print(file = bashFile)
#
#  # Generate a vcf of X de novo variants
#  print("  # X de novo variants", file = bashFile)
#  print("  echo \"Generating vcf of x de novo variants...\"", file = bashFile)
#  print("  X_DENOVO_VCF=", moiFiles["xdenovoVcf"]["file"], sep = "", file = bashFile)
#  print("  bcftools view -i 'INFO/x_denovo=\"", proband, "\"' -O z -o $X_DENOVO_VCF $FILTEREDVCF", sep = "", file = bashFile)
#  print(file = bashFile)
#
#  # Generate a vcf of de novo variants
#  print("  # De novo variants", file = bashFile)
#  print("  echo \"Generating vcf of recessive variants...\"", file = bashFile)
#  print("  RECESSIVE_VCF=", moiFiles["recessiveVcf"]["file"], sep = "", file = bashFile)
#  print("  bcftools view -i 'INFO/recessive=\"", proband, "\"' -O z -o $RECESSIVE_VCF $FILTEREDVCF", sep = "", file = bashFile)
#  print(file = bashFile)
#
#  # Generate a vcf of comp het variants
#  print("  # Compound heterozygous variants", file = bashFile)
#  print("  echo \"Generating vcf of compound heterozygous variants...\"", file = bashFile)
#  print("  COMPHET_VCF=", moiFiles["comphetVcf"]["file"], sep = "", file = bashFile)
#  print("  bcftools view -i 'INFO/slivar_comphet=\"", proband, "\"' -O z -o $COMPHET_VCF $FILTEREDVCF", sep = "", file = bashFile)
#  print(file = bashFile)
#
## Generate the tsv files to pass annotations to Mosaic
#def generateTsv(args, bashFile):
#  global resourceInfo
#  global mosaicInfo
#  privateResources = []
#
#  # Loop over all the annotations to pass to Mosaic. This is only the annotations in the Mosaic resource file
#  # which may contain annotations not in the resources file. This can happen if the annotation is not generated
#  # from a resource file, e.g. HPO terms.
#  for resource in resourceInfo["resources"]:
#  #for resource in mosaicInfo["resources"]:
#
#    # Check if this resource is to be uploaded to Mosaic, and if it can be (e.g. does it appear
#    # in the Mosaic json
#    upload    = resourceInfo["resources"][resource]["upload"]
#    canUpload = True if resource in mosaicInfo["resources"] else False
#
#    # If the resource is listed as to be uploaded to Mosaic, but there are no instructions on how
#    # to, fail
#    if upload and not canUpload: fail("Resource '" + str(resource) + "' is marked as to be uploaded to Mosaic, but it is not included in the Mosaic json")
#
#    # If the resource is not listed as to be uploaded to Mosaic, but it can be, provide a warning.
#    elif not upload and canUpload:
#      warningTitle = "Resource omitted"
#      description  = "The following resources were listed as to be uploaded to Mosaic in the resources json, "
#      description += "but no instructions are provided in the Mosaic json, so they will not be uploaded"
#      if warningTitle not in warnings: warnings[warningTitle] = {"description": description, "messages": [resource]}
#      else: warnings[warningTitle]["messages"].append(resource)
#
#    # For resources to upload, check if they are public or private Mosaic annotations. All private
#    # annotations should be included in the same tsv, but public annotation will get their own tsv.
#    elif upload and canUpload:
#      annotation_type = mosaicInfo["resources"][resource]["annotation_type"]
#      if annotation_type == "private": privateResources.append(resource)
#      else: buildPublicTsv(args, resource, bashFile)
#
#  # Build the tsv for private annotations
#  if len(privateResources) > 0:
#
#    # All of the private annotations need to be created in the project
#    temp = 1
#    annotationUids = createAnnotations(args, privateResources)
#    buildPrivateAnnotationTsv(args, privateResources, annotationUids, bashFile)
#
## Build the tsv file
#def buildPublicTsv(args, resource, bashFile):
#  global mosaicInfo
#  global tsvFiles
#
#  # First create the header line
#  header = "CHROM\\tSTART\\tEND\\tREF\\tALT"
#  for annotation in sorted(mosaicInfo["resources"][resource]["annotations"]):
#    header += "\\t" + mosaicInfo["resources"][resource]["annotations"][annotation]["uid"]
#
#  # Print the header
#  print("  # Public annotation: ", resource, sep = "", file = bashFile)
#  outputFile         = resource + "_calypso_annotations_mosaic.tsv"
#  tsvFiles[resource] = outputFile
#  tempFile           = outputFile + ".tmp"
#  print("  echo -e \"", header, "\" > ", outputFile, sep = "", file = bashFile)
#
#  # Include all annotations for all resources
#  annotations = "%CHROM\\t%POS\\t%END\\t%REF\\t%ALT"
#
#  delimiter     = mosaicInfo["resources"][resource]["delimiter"]
#  isPostprocess = mosaicInfo["resources"][resource]["postprocess"]
#
#  # Loop over the annotations and add to the command line
#  for annotation in sorted(mosaicInfo["resources"][resource]["annotations"]):
#    if delimiter: annotations += "\\t%INFO/" + mosaicInfo["resources"][resource]["info"]
#    else: annotations += "\\t%INFO/" + annotation
#
#  # Build the query command
#  print("  bcftools query -f '", annotations, "\\n' \\", sep = "", file = bashFile)
#  print("  $FILTEREDVCF >> ", outputFile, sep = "", file = bashFile)
#  print(file = bashFile)
#
#  # If the resources need post processing (e.g. with additional scripts), do not modify the tsv
#  if isPostprocess:
#    print("  # Postprocess the annotation tsv", file = bashFile)
#    postCommand = mosaicInfo["resources"][resource]["post_command"]
#    command = (postCommand["precommand"] + " ") if postCommand["precommand"] else ""
#    command += os.path.dirname( __file__) + "/" + postCommand["tool"]
#    for argument in postCommand["args"]:
#      argValue = postCommand["args"][argument]
#
#      # Process the different arguments for standard values
#      if argValue == "TSV": argValue = outputFile
#      elif argValue == "data_directory": argValue = args.data_directory
#      elif argValue == "reference": argValue = args.reference
#      command += " " + argument + " " + argValue
#    print("  ", command, sep = "", file = bashFile)
#    print(file = bashFile)
#  else:
#
#    # If the delimiter has been set, the tsv contains the full annotation (e.g. A|B|C), for each
#    # annotation. This needs to be modified to break the annotation up.
#    print("  # Update the tsv file to Mosaic specifications", file = bashFile)
##    print("  awk '{FS=\"\\t\"; OFS=\"\\t\"} {if ($1 != \"CHROM\") {if ($1 ~ \"chr\") {$1=substr($1, 4)}; $3 = $3+1; " + awkCommand + "for(i=6; i<=NF; i++) {$i=($i==\".\" ? \"\":$i)}}; print $0}' \\", file = bashFile)
#    print("  awk '{FS=\"\\t\"; OFS=\"\\t\"} {if ($1 != \"CHROM\") {if ($1 ~ \"chr\") {$1=substr($1, 4)}; $3 = $3+1; for(i=6; i<=NF; i++) {$i=($i==\".\" ? \"\":$i)}}; print $0}' \\", file = bashFile)
#    print(" ", outputFile, ">", tempFile, file = bashFile)
#    print("  mv", tempFile, outputFile, file = bashFile)
#    print(file = bashFile)
#
## Creat new annotations
#def createAnnotations(args, privateResources):
#  global mosaicConfig
#  global mosaicInfo
#  global samples
#  global sampleOrder
#  global createdAnnotations
#  annotationUids = {}
#
#  # Get all the annotations in the project
#  data = json.loads(os.popen(api_va.getVariantAnnotations(mosaicConfig, args.project_id)).read())
#
#  # Store the annotations in the project by name
#  projectAnnotations = {}
#  for annotation in data: projectAnnotations[str(annotation["name"])] = annotation
#
#  # Loop over the private resources
#  for resource in sorted(privateResources):
#    annotationUids[resource] = {}
#
#    # Loop over all the annotations for this resource
#    for annotation in mosaicInfo["resources"][resource]["annotations"]:
#
#      # If this is a genotype annotation, create an attribute for each sample
#      isGenotype = mosaicInfo["resources"][resource]["annotations"][annotation]["isGenotype"]
#      annotationUids[resource][annotation] = {"isGenotype": isGenotype, "uids": []}
#      if isGenotype:
#
#        # Loop over all samples in the project and get the annotations for each sample
#        for sample in sampleOrder:
#          annotationName = str(annotation) + " " + str(samples[sample]["relationship"])
#          annotationUid  = checkPrivateAnnotation(args, resource, annotation, projectAnnotations, annotationName, annotationUids)
#
#      # Other types of private annotations are not yet handled
#      else:
#        annotationUid = checkPrivateAnnotation(args, resource, annotation, projectAnnotations, annotation, annotationUids)
#
#  return annotationUids
#
## Check if a private annotation exists and create if it doesn't
#def checkPrivateAnnotation(args, resource, annotation, projectAnnotations, annotationName, annotationUids):
#  global mosaicConfig
#  global mosaicInfo
#  global createdAnnotations
#
#  # Check if an annotation with this name already exists
#  annotationUid = projectAnnotations[annotationName]["uid"] if annotationName in projectAnnotations else False
#
#  # If the annotation does not exist, create it
#  if not annotationUid:
#    valueType = mosaicInfo["resources"][resource]["annotations"][annotation]["type"]
#    data      = json.loads(os.popen(api_va.postCreateVariantAnnotations(mosaicConfig, annotationName, valueType, "private", args.project_id)).read())
#
#    # Get the uid of the newly created annotation
#    annotationUid = data["uid"]
#
#  # Store the annotation uids
#  annotationUids[resource][annotation]["uids"].append(annotationUid)
#  createdAnnotations[annotationName] = annotationUid
#
#  return annotationUids
#
## Build the tsv containing the private annotations
#def buildPrivateAnnotationTsv(args, resources, annotationUids, bashFile):
#  global mosaicConfig
#  global tsvFiles
#  global sampleOrder
#
#  # First create the header line
#  header = "CHROM\\tSTART\\tEND\\tREF\\tALT"
#  for resource in sorted(resources):
#    for annotation in sorted(annotationUids[resource]):
#      for uid in annotationUids[resource][annotation]["uids"]: header += "\\t" + uid
#
#  # Print the header
#  print("  # Private annotation: ", resource, sep = "", file = bashFile)
#  outputFile         = resource + "_calypso_annotations_mosaic.tsv"
#  tsvFiles[resource] = outputFile
#  tempFile           = outputFile + ".tmp"
#  print("  echo -e \"", header, "\" > ", outputFile, sep = "", file = bashFile)
#
#  # Include all annotations for all resources
#  annotations = "%CHROM\\t%POS\\t%END\\t%REF\\t%ALT"
#
#  delimiter     = mosaicInfo["resources"][resource]["delimiter"]
#  isPostprocess = mosaicInfo["resources"][resource]["postprocess"]
#
#  # Loop over the annotations and add to the command line
#  for annotation in sorted(mosaicInfo["resources"][resource]["annotations"]):
#
#    # Determine if this is a genotype annotation, loop over the samples in the correct order, and add the annotation
#    if mosaicInfo["resources"][resource]["annotations"][annotation]["isGenotype"]: annotations += "[\\t%GQ]"
#    elif delimiter: annotations += "\\t%INFO/" + mosaicInfo["resources"][resource]["info"]
#    else: annotations += "\\t%INFO/" + annotation
#
#  # Build the query command
#  print("  bcftools query -f '", annotations, "\\n' \\", sep = "", file = bashFile)
#  print("  $FILTEREDVCF >> ", outputFile, sep = "", file = bashFile)
#  print(file = bashFile)
#
#  # If the resources need post processing (e.g. with additional scripts), do not modify the tsv
#  if isPostprocess:
#    print("  # Postprocess the annotation tsv", file = bashFile)
#    postCommand = mosaicInfo["resources"][resource]["post_command"]
#    command = (postCommand["precommand"] + " ") if postCommand["precommand"] else ""
#    command += os.path.dirname( __file__) + "/" + postCommand["tool"]
#    for argument in postCommand["args"]:
#      argValue = postCommand["args"][argument]
#
#      # Process the different arguments for standard values
#      if argValue == "TSV": argValue = outputFile
#      elif argValue == "data_directory": argValue = args.data_directory
#      elif argValue == "reference": argValue = args.reference
#      command += " " + argument + " " + argValue
#    print("  ", command, sep = "", file = bashFile)
#    print(file = bashFile)
#  else:
#
#    # If the delimiter has been set, the tsv contains the full annotation (e.g. A|B|C), for each
#    # annotation. This needs to be modified to break the annotation up.
#    print("  # Update the tsv file to Mosaic specifications", file = bashFile)
##    print("  awk '{FS=\"\\t\"; OFS=\"\\t\"} {if ($1 != \"CHROM\") {if ($1 ~ \"chr\") {$1=substr($1, 4)}; $3 = $3+1; " + awkCommand + "for(i=6; i<=NF; i++) {$i=($i==\".\" ? \"\":$i)}}; print $0}' \\", file = bashFile)
#    print("  awk '{FS=\"\\t\"; OFS=\"\\t\"} {if ($1 != \"CHROM\") {if ($1 ~ \"chr\") {$1=substr($1, 4)}; $3 = $3+1; for(i=6; i<=NF; i++) {$i=($i==\".\" ? \"\":$i)}}; print $0}' \\", file = bashFile)
#    print(" ", outputFile, ">", tempFile, file = bashFile)
#    print("  mv", tempFile, outputFile, file = bashFile)
#    print(file = bashFile)
#
## If HPO ids were supplied, process them
#def processHpo(args, bashFile):
#  global resourceInfo
#  global mosaicInfo
#  global rootPath
#
#  hpoTermsUid    = mosaicInfo['resources']['hpo']['annotations']['HPO Terms']['uid']
#  hpoOverlapsUid = mosaicInfo['resources']['hpo']['annotations']['HPO Overlaps']['uid']
#  print("# Processing HPO terms...", file = bashFile)
#  print("  echo -e \"Processing HPO terms to generate additional HPO annotations...\"", file = bashFile)
#  print("  bcftools query -f '%CHROM\\t%POS\\t%END\\t%REF\\t%ALT\\t%INFO/BCSQ\\n' \\", file = bashFile)
#  print("  $FILTEREDVCF \\", file = bashFile)
#  print("  | python ", rootPath, "/calypso_hpo.py \\", sep = "", file = bashFile)
#  print("  -p ", args.project_id, " \\", sep = "", file = bashFile)
#  print("  -c \"", args.config, "\" \\", sep = "", file = bashFile)
#  print("  -o \"", resourceInfo["resources"]['hpo']['file'], "\" \\", sep = "", file = bashFile)
#  print("  -r \"", args.hpo, "\" \\", sep = "", file = bashFile)
#  print("  -e \"", hpoTermsUid, "\" \\", sep = "", file = bashFile)
#  print("  -l \"", hpoOverlapsUid, "\"", sep = "", file = bashFile)
#  print(file = bashFile)
#
## Check the public annotations exist, and get their ids
#def getPublicAnnotation(args, resource):
#  global mosaicInfo
#  global warnings
#  page = 1
#  uids = {}
#
#  # Find the annotation id for the public annotation. Get the first page of annotations and check
#  # if there are additional pages required
#  for annotation in mosaicInfo["resources"][resource]["annotations"]:
#    uids[mosaicInfo["resources"][resource]["annotations"][annotation]["uid"]] = False
#
#  # Get the first page of public attributes, and how many pages there are
#  noPages          = getAnnotationsPages(args)
#  isComplete, uids = getAnnotations(args, page, uids)
#
#  # If some ids are not found, go to the next page. If there are no more pages, throw a warning
#  while page <= noPages:
#    isComplete, uids = getAnnotations(args, page, uids)
#    page += 1
#    if isComplete: break
#
#  # If the annotations weren't found, throw a warning
#  if not isComplete:
#    warningTitle = "Annotation not found"
#    description  = "Public annotations need to be imported into the Mosaic project, but these annotations weren't found "
#    description += "and so cannot be imported"
#    if warningTitle not in warnings: warnings[warningTitle] = {"description": description, "messages": []}
#    for annotation in uids:
#      if not uids[annotation]: warnings[warningTitle]["messages"].append("  " + annotation)
#
#  # If the annotations can be imported, import them
#  else:
#    for annotation in uids: importAnnotation(args, annotation, uids[annotation])
#
#  # Return the isComplete variable
#  return isComplete
#
## Determine how many pages of public attributes there are
#def getAnnotationsPages(args):
#  global mosaicConfig
#
#  # Get a page of public annotations from Mosaic (100 annotations per page)
#  data = json.loads(os.popen(api_va.getVariantAnnotationsImport(mosaicConfig, args.project_id, 100, 1)).read())
#  noPages = math.ceil(int(data["count"]) / int(100))
#
#  # Return the number of pages of annotations
#  return noPages
#
## Get a page of public attributes
#def getAnnotations(args, page, uids):
#  global mosaicConfig
#
#  # Get a page of public annotations from Mosaic (100 annotations per page)
#  data = json.loads(os.popen(api_va.getVariantAnnotationsImport(mosaicConfig, args.project_id, 100, page)).read())
#  for annotation in data["data"]:
#    if annotation["uid"] in uids.keys(): uids[annotation["uid"]] = annotation["id"]
#
#  # Check if all required attributes have been found and the id stored
#  isComplete = True
#  for annotation in uids:
#    if not uids[annotation]: isComplete = False
#
#  # Return the Mosaic uids and whether they were all found
#  return isComplete, uids
#
## Import a public annotation into the Mosaic project
#def importAnnotation(args, annotation, uid):
#  global mosaicConfig
#
#  # Create the command to import the annotation
#  try: data = json.loads(os.popen(api_va.postImportVariantAnnotations(mosaicConfig, uid, args.project_id)).read())
#  except: fail("Unable to import variant annotation: " + uid)
#
## Output summary file
#def calypsoSummary(args, finalVcf):
#  global resourceInfo
#  global workingDir
#  global version
#  global date
#
#  # Open an output summary file
#  summaryFileName  = workingDir + "calypso_" + date + ".txt"
#  try: summaryFile = open(summaryFileName, "w")
#  except: fail("Failed to open summary file")
#
#  # Write relevant information to file
#  print("### Output from Calypso pipeline ###", file = summaryFile)
#  print(file = summaryFile)
#  print("Calypso pipeline version: ", version, sep = "", file = summaryFile)
#  print("Calypso resource version: ", resourceInfo['version'], sep = "", file = summaryFile)
#  print("Reference:                ", args.reference, sep = "", file = summaryFile)
#  print("Created on:               ", date, sep = "", file = summaryFile)
#  print("Generated VCF file:       ", finalVcf, sep = "", file = summaryFile)
#  print(file = summaryFile)
#
#  # Loop over all the used resources and output their versions
#  for resource in resourceInfo["resources"]:
#    print(resource, file = summaryFile)
#    print("  version: ", resourceInfo["resources"][resource]["version"], sep = "", file = summaryFile)
#    print("  file:    ", resourceInfo["resources"][resource]["file"], sep = "", file = summaryFile)
#    print(file = summaryFile)
#
#  # Close the file
#  summaryFile.close()
#
## Update the Calypso history value
#def updateHistory(args, valueRecords):
#  global version
#  global date
#
#  # Get the history string for this project
#  for record in valueRecords:
#    value = record["value"] if str(record["project_id"]) == str(args.project_id) else False
#
#  # If there was no history value, return todays date and version
#  if not value: return "v" + str(version) + ":" + date
#
#  # Get the most recent update to the history
#  mostRecent  = value.split(",")[-1] if "," in value else value
#  lastVersion = mostRecent.split(":")[0]
#  lastDate    = mostRecent.split(":")[1]
#
#  # If the date and version at execution are the same as at the last execution, do not update
#  if str(lastVersion) == str(version) and str(lastDate) == date: return value
#  else: return str(value) + "," + str(version) + ":" + date

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)

# Initialise global variables

# Pipeline version
version = "1.0.0"
date    = str(date.today())

# The working directory where all created files are kept
workingDir = os.getcwd() + "/calypso_v" + version + "r"

## Store warnings to be output at the end
#warnings = {}
#
## Store info on allowed values
#isFamily          = True
allowedFamily     = ['singleton', 'duo', 'trio', 'quad', 'quintet']
allowedReferences = ['37', '38']
#
## Store the mode of inheritance file names
#moiFiles = {}

# Store information related to Mosaic
mosaicConfig = {}

## Store the tsv files that upload annotations to Mosaic
#tsvFiles = {}
#
## Store the path to the Calypso directory
#rootPath = os.path.dirname( __file__)
#
## Store all the annotations in a project
#projectAnnotations = {}
#
## Store annotations created for this project
#createdAnnotations = {}
#
## Created files
#tomlFilename = ""

if __name__ == "__main__":
  main()
