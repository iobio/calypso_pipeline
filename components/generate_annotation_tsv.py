#!/usr/bin/python

from __future__ import print_function
from datetime import date
from os.path import exists
from sys import path

import argparse
import os
import sys

# Add the path of the common functions and import them
path.append(os.path.dirname(__file__) + '/components')
import tools_bcftools as bcftools
import calypso_mosaic_resources as mosr

def main():
  global allowedClasses

  # Parse the command line
  args = parseCommandLine()

  # Define the executable bcftools command
  if args.tools_directory:
    if not args.tools_directory.endswith('/'): args.tools_directory = str(args.tools_directory) + '/'
    bcftoolsExe = args.tools_directory + 'bcftools/bcftools'
  else: bcftoolsExe = 'bcftools'

  # Read the mosaicJson file to get information on how to process different annotations
  mosaicInfo = mosr.readMosaicJson(args.mosaic_json, args.reference)

  # Check that the provided annotation class is valid
  processClass = mosaicInfo['resources'][args.resource]['class']
  if processClass in allowedClasses:

    # Open an output tsv file to write annotations to
    outputFile = open(args.output_tsv, 'w')
  
    # Write the header line to the tsv file
    print('CHROM\tSTART\tEND\tREF\tALT\t', '\t'.join(args.uids.replace(' ', '').split(',')), sep = '', file = outputFile)
  
    # Loop over the vcf and process according to the annotation class
    if processClass == 'A': processClassA(args.input_vcf, args.tags.replace(' ', '').split(','), outputFile)
    elif processClass == 'B': processClassB(args.input_vcf, args.tags.replace(' ', '').split(','), outputFile)
    elif processClass == 'C': processClassC(mosaicInfo, args.resource, args.input_vcf, args.tags.replace(' ', '').split(','), outputFile)
    elif processClass == 'OMIM': processClassOMIM(args.input_vcf, args.tags.replace(' ', '').split(','), outputFile)
    elif processClass == 'compound': processClassCompound(args.resource, mosaicInfo['resources'][args.resource], args.input_vcf, args.tags.replace(' ', '').split(','), outputFile)
    elif processClass == 'spliceai': processClassSpliceAI(mosaicInfo['resources'][args.resource], args.input_vcf, args.tags.split(','), outputFile)
    elif processClass == 'clinvar': processClassClinvar(resourceInfo, args.input_vcf, args.tags.replace(' ', '').split(','), outputFile)
  
    # Close the output tsv file
    outputFile.close()

  # Write a warning if the annotation class is not recognised and do not process the vcf file
  else: print('Unable to process annotations for resource: ' + args.resource, sep = '', file = sys.stderr)

# Input options
def parseCommandLine():
  global version
  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  parser.add_argument('--config', '-c', required = True, metavar = "string", help = "The config file for Mosaic")
  parser.add_argument('--input_vcf', '-i', required = True, metavar = 'string', help = 'The input vcf file to annotate')
  parser.add_argument('--output_tsv', '-o', required = True, metavar = 'string', help = 'The output tsv file')
  parser.add_argument('--resource', '-e', required = True, metavar = 'string', help = 'The name of the resource (used for the output tsv name')
  parser.add_argument('--reference', '-r', required = True, metavar = 'string', help = 'The genome reference file used')
  parser.add_argument('--tags', '-g', required = True, metavar = 'string', help = 'A comma separated list of VCF INFO tags to be extracted')
  parser.add_argument('--uids', '-d', required = True, metavar = 'string', help = 'A comma separated list of uids for the resource annotations')
  parser.add_argument('--tools_directory', '-s', required = False, metavar = 'string', help = 'The path to the directory where the tools live')
  parser.add_argument('--mosaic_json', '-m', required = True, metavar = 'string', help = 'The json file describing the Mosaic parameters')

  # Optional mosaic arguments
  parser.add_argument('--token', '-t', required = False, metavar = "string", help = "The Mosaic authorization token")
  parser.add_argument('--url', '-u', required = False, metavar = "string", help = "The base url for Mosaic")
  parser.add_argument('--attributes_project', '-a', required = False, metavar = "integer", help = "The Mosaic project id that contains public attributes")

  # Version
  parser.add_argument('--version', '-v', action="version", version='Calypso annotation pipeline version: ' + str(version))

  return parser.parse_args()

# Process class A annotations. This is for annotations that are floats and if multiple values occur, only the maximum value should be
# uploaded to Mosaic. This includes the following annoatations:
#   - CCR
#   - EVE
#   - gnomAD
#   - MutScore
#   - pLI
#   - REVEL
def processClassA(vcf, tags, outputFile):

  # Loop over all records in the vcf file
  for record in os.popen(bcftools.query(bcftoolsExe, vcf, tags)).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')
  
    # Check that the annotation has a value. If not, skip this record. The only records that will be output are those with values
    hasValue = False
    for i in range(5, len(fields)):
      if fields[5] != '.':
        hasValue = True
        break
    if hasValue:
  
      # Update the chromosome and position
      fields[0], fields[2] = updateCoords(fields[0], fields[2])
  
      # If there are multiple values, only output the largest
      for i in range(5, len(fields)):
        if ',' in fields[i]:
          outputValue = 0.
          for value in fields[i].split(','):
            if float(value) > float(outputValue): outputValue = value
          fields[i] = str(outputValue)
  
      # Check that the value is a float
      try: typeTest = float(fields[5])
      except: fail('Invalid value for annotation: ' + record.rstrip())
  
      # Make sure the value falls between 1E-37 and 1E+37
      if float(fields[5]) <= 1e-37: fields[5] = "1e-37"
      elif float(fields[5]) >= 1e37: fields[5] = "1e37"
  
      # Build the output record from the updated fields
      print('\t'.join(fields), file = outputFile)

# Process class B annotations. This is for annotations that are strings that do not undergo any modifications
#   - dbSNP
def processClassB(vcf, tags, outputFile):

  # Loop over all records in the vcf file
  for record in os.popen(bcftools.query(bcftoolsExe, vcf, tags)).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')
  
    # Check that the variant has a value. If not, skip this record. The only records that will be output are those with values
    hasValue = False
    for i in range(5, len(fields)):
      if fields[5] != '.':
        hasValue = True
        break
    if hasValue:
  
      # Update the chromosome and position
      fields[0], fields[2] = updateCoords(fields[0], fields[2])
  
      # Build the output record from the updated fields
      print('\t'.join(fields), file = outputFile)

# Process class C annotations. This is for annotations that are found in the genotype (FORMAT) strings in the vcf
#   - GQ
def processClassC(mosaicInfo, resource, vcf, tags, outputFile):

  # Loop over all records in the vcf file
  for record in os.popen(bcftools.queryFormat(bcftoolsExe, vcf, resource)).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')
  
    # If all genotype fields are '.', do not output the line. If some have values and others are '.', replace the '.' with ' '
    hasValue = False
    for i, annotation in enumerate(fields[5:]):
      if annotation == '.': fields[i + 5] = ''
      else: hasValue = True

    # If any of the fields have values, update and output
    if hasValue:

      # Update the chromosome and position
      fields[0], fields[2] = updateCoords(fields[0], fields[2])

      # Build the output record from the updated fields
      print('\t'.join(fields), file = outputFile)

# Process compound annotations
def processClassCompound(resource, resourceInfo, vcf, tags, outputFile):

  # Check that the delimeter is provided to determine how to split up the compound annotation. If it is not set, then provide a
  # warning and do not preceed with this annotation
  try: delimeter = resourceInfo['delimeter']
  except:
    print('The delimeter field is not provided for resource ', resource, ' and so its annotation cannot be processed.', sep = '')
    return

  # Loop over all records in the vcf file
  for record in os.popen(bcftools.splitvepQuery(bcftoolsExe, vcf, [resourceInfo['info_field']])).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')
  
    # Check that the variant has a value. If not, skip this record. The only records that will be output are those with values
    if fields[5] != '.':
  
      # Update the chromosome and position
      fields[0], fields[2] = updateCoords(fields[0], fields[2])
      annotations = fields.pop().split(delimeter)
      hasValue = False
      for tag in tags:

        # If the annotation is longer than 255 characters, it cannot be uploaded to Mosaic, so trim it
        value = annotations[resourceInfo['annotations'][tag]['position'] - 1]
        if len(value) > 255: value = value[0:252] + '...'
        fields.append(value)
        if value: hasValue = True
  
      # Build the output record from the updated fields
      if hasValue: print('\t'.join(fields), file = outputFile)

# Process OMIM annotations
def processClassOMIM(vcf, tags, outputFile):

  # Loop over all records in the vcf file
  for record in os.popen(bcftools.query(bcftoolsExe, vcf, tags)).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')
  
    # Check that the variant has a value. If not, skip this record. The only records that will be output are those with values
    hasValue = False
    for i in range(5, len(fields)):
      if fields[5] != '.':
        hasValue = True
        break
    if hasValue:
  
      # Update the chromosome and position
      fields[0], fields[2] = updateCoords(fields[0], fields[2])
  
      # If the field begins with a comma, or has two consecutive commas, there is no value. For example, a variant may be associated
      # with three MIM ids, but only the third has an associated inheritance. The inheritance will appear as ',,AR', for example. In
      # this case, the leading ',' should be replaced with a leading 'None', and the consecutive commas should be replaced with
      # ',None,' yielding a value of None,None,AR for upload to Mosaic
      for i in range(5, len(fields)):
        if fields[i].startswith(','): fields[i] = 'None' + fields[i]
        fields[i] = fields[i].replace(',,', ',None,')
  
      # Build the output record from the updated fields
      print('\t'.join(fields), file = outputFile)

# Process SpliceAI annotations
def processClassSpliceAI(resourceInfo, vcf, tags, outputFile):

  # Loop over all records in the vcf file
  for record in os.popen(bcftools.query(bcftoolsExe, vcf, [resourceInfo['info_field']])).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')
  
    # Check that the variant has a value. If not, skip this record. The only records that will be output are those with values
    if fields[5] != '.':
  
      # Update the chromosome and position
      fields[0], fields[2] = updateCoords(fields[0], fields[2])
      annotations = fields.pop()

      # If there are multiple sets of SpliceAI values, output the one with the largest score
      spliceAI = annotations.split(',') if ',' in annotations else [annotations]
      maxScore = -1.

      # Loop over the available sets of scores
      finalSet = []
      for annotation in spliceAI:
        annotations = annotation.split('|')
        temp        = []
        isMax       = False
        for tag in tags:
          value = annotations[resourceInfo['annotations'][tag]['position'] - 1]
  
          # If this is the SpliceAI Max Score, take the max value from the supplied positions
          if tag == 'SpliceAI Max Score':
            values = []
            for i in resourceInfo['annotations'][tag]['positions']: values.append(annotations[i - 1])
            value = max(values)
            try: typeTest = float(value)
            except: fail('Invalid value for annotation: ' + record.rstrip())
            temp.append(value)

            # If the max SpliceAI is greater than what has been seen to date, this is the set of scores to be output
            if float(value) > float(maxScore): isMax = True
          else: 
            value = annotations[resourceInfo['annotations'][tag]['position'] - 1]
            try: typeTest = float(value)
            except: fail('Invalid value for annotation: ' + record.rstrip())
            temp.append(value)
 
        # If this set of scores had the largest max score to date, store this set
        if isMax: finalSet = temp
    
      # Build the output record from the updated fields
      print('\t'.join(fields + finalSet), file = outputFile)

# Process ClinVar annotations
def processClassClinvar(resourceInfo, vcf, tags, outputFile):
  skippedValues    = {}
  noRecords        = 0
  noOutputRecords  = 0
  noSkippedRecords = 0
  noDotRecords     = 0

  # Define the allowed ClinVar significance values
  allowedClinVar = {}
  allowedClinVar['Benign'] = 'Benign'
  allowedClinVar['Benign/Likely_benign'] = 'Benign/Likely_benign'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity'] = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Likely_benign'] = 'Likely_benign'
  allowedClinVar['Pathogenic'] = 'Pathogenic'
  allowedClinVar['Pathogenic/Likely_pathogenic'] = 'Pathogenic/Likely_pathogenic'
  allowedClinVar['Likely_pathogenic'] = 'Likely_pathogenic'
  allowedClinVar['Uncertain_significance'] = 'Uncertain_significance'

  # Modify some of the obsevered values
  allowedClinVar['Pathogenic,_risk_factor']                            = 'Pathogenic'
  allowedClinVar['Pathogenic|risk_factor']                             = 'Pathogenic'
  allowedClinVar['Pathogenic|_risk_factor']                            = 'Pathogenic'
  allowedClinVar['Pathogenic,_other']                                  = 'Pathogenic'
  allowedClinVar['Pathogenic|other']                                   = 'Pathogenic'
  allowedClinVar['Pathogenic|_other']                                  = 'Pathogenic'
  allowedClinVar['Pathogenic,_Affects']                                = 'Pathogenic'
  allowedClinVar['Pathogenic|Affects']                                 = 'Pathogenic'
  allowedClinVar['Pathogenic|_Affects']                                = 'Pathogenic'
  allowedClinVar['Pathogenic,_protective']                             = 'Pathogenic'
  allowedClinVar['Pathogenic|protective']                              = 'Pathogenic'
  allowedClinVar['Pathogenic|_protective']                             = 'Pathogenic'
  allowedClinVar['Pathogenic,_drug_response']                          = 'Pathogenic'
  allowedClinVar['Pathogenic|drug_response']                           = 'Pathogenic'
  allowedClinVar['Pathogenic|_drug_response']                          = 'Pathogenic'
  allowedClinVar['Pathogenic,_association']                            = 'Pathogenic'
  allowedClinVar['Pathogenic|association']                             = 'Pathogenic'
  allowedClinVar['Pathogenic|_association']                            = 'Pathogenic'
  allowedClinVar['Pathogenic,_association,_protective']                = 'Pathogenic'
  allowedClinVar['Pathogenic|association|protective']                  = 'Pathogenic'
  allowedClinVar['Pathogenic|_association|_protective']                = 'Pathogenic'
  allowedClinVar['Pathogenic,_drug_response,_other']                   = 'Pathogenic'
  allowedClinVar['Pathogenic|drug_response|other']                     = 'Pathogenic'
  allowedClinVar['Pathogenic|_drug_response|_other']                   = 'Pathogenic'
  allowedClinVar['Pathogenic,_confers_sensitivity']                    = 'Pathogenic'
  allowedClinVar['Pathogenic|confers_sensitivity']                     = 'Pathogenic'
  allowedClinVar['Pathogenic|_confers_sensitivity']                    = 'Pathogenic'
  allowedClinVar['Pathogenic,_drug_response,_protective,_risk_factor'] = 'Pathogenic'
  allowedClinVar['Pathogenic/Likely_risk_allele']                      = 'Pathogenic'
  allowedClinVar['Pathogenic/Pathogenic,_low_penetrance']              = 'Pathogenic'

  allowedClinVar['Pathogenic/Likely_pathogenic,_other']                           = 'Pathogenic/Likely_pathogenic'
  allowedClinVar['Pathogenic/Likely_pathogenic|other']                            = 'Pathogenic/Likely_pathogenic'
  allowedClinVar['Pathogenic/Likely_pathogenic|_other']                           = 'Pathogenic/Likely_pathogenic'
  allowedClinVar['Pathogenic/Likely_pathogenic,_risk_factor']                     = 'Pathogenic/Likely_pathogenic'
  allowedClinVar['Pathogenic/Likely_pathogenic|risk_factor']                      = 'Pathogenic/Likely_pathogenic'
  allowedClinVar['Pathogenic/Likely_pathogenic|_risk_factor']                     = 'Pathogenic/Likely_pathogenic'
  allowedClinVar['Pathogenic/Likely_pathogenic,_drug_response']                   = 'Pathogenic/Likely_pathogenic'
  allowedClinVar['Pathogenic/Likely_pathogenic|drug_response']                    = 'Pathogenic/Likely_pathogenic'
  allowedClinVar['Pathogenic/Likely_pathogenic|_drug_response']                   = 'Pathogenic/Likely_pathogenic'
  allowedClinVar['Pathogenic/Pathogenic,_low_penetrance|risk_factor']             = 'Pathogenic/Likely_pathogenic'
  allowedClinVar['Pathogenic/Likely_pathogenic/Likely_risk_allele']               = 'Pathogenic/Likely_pathogenic'
  allowedClinVar['Pathogenic/Likely_pathogenic,_other,_risk_factor']              = 'Pathogenic/Likely_pathogenic'
  allowedClinVar['Pathogenic/Likely_pathogenic,_association']                     = 'Pathogenic/Likely_pathogenic'
  allowedClinVar['Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance']       = 'Pathogenic/Likely_pathogenic'
  allowedClinVar['Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance|other'] = 'Pathogenic/Likely_pathogenic'
  allowedClinVar['Likely_pathogenic/Pathogenic,_low_penetrance']                  = 'Pathogenic/Likely_pathogenic'

  allowedClinVar['Likely_pathogenic,_other']              = 'Likely_pathogenic'
  allowedClinVar['Likely_pathogenic|other']               = 'Likely_pathogenic'
  allowedClinVar['Likely_pathogenic|_other']              = 'Likely_pathogenic'
  allowedClinVar['Likely_pathogenic,_risk_factor']        = 'Likely_pathogenic'
  allowedClinVar['Likely_pathogenic|risk_factor']         = 'Likely_pathogenic'
  allowedClinVar['Likely_pathogenic|_risk_factor']        = 'Likely_pathogenic'
  allowedClinVar['Likely_pathogenic,_association']        = 'Likely_pathogenic'
  allowedClinVar['Likely_pathogenic|association']         = 'Likely_pathogenic'
  allowedClinVar['Likely_pathogenic|_association']        = 'Likely_pathogenic'
  allowedClinVar['Likely_pathogenic,_drug_response']      = 'Likely_pathogenic'
  allowedClinVar['Likely_pathogenic|drug_response']       = 'Likely_pathogenic'
  allowedClinVar['Likely_pathogenic|_drug_response']      = 'Likely_pathogenic'
  allowedClinVar['Likely_pathogenic,_other,_risk_factor'] = 'Likely_pathogenic'
  allowedClinVar['Likely_pathogenic,_Affects']            = 'Likely_pathogenic'
  allowedClinVar['Likely_pathogenic|Affects']             = 'Likely_pathogenic'
  allowedClinVar['Likely_pathogenic|_Affects']            = 'Likely_pathogenic'
  allowedClinVar['Likely_pathogenic/Likely_risk_allele']  = 'Likely_pathogenic'
  allowedClinVar['Likely_pathogenic,_low_penetrance']     = 'Likely_pathogenic'

  allowedClinVar['Uncertain_significance,_risk_factor']          = 'Uncertain_significance'
  allowedClinVar['Uncertain_significance|risk_factor']           = 'Uncertain_significance'
  allowedClinVar['Uncertain_significance|_risk_factor']          = 'Uncertain_significance'
  allowedClinVar['Uncertain_significance,_other']                = 'Uncertain_significance'
  allowedClinVar['Uncertain_significance|other']                 = 'Uncertain_significance'
  allowedClinVar['Uncertain_significance|_other']                = 'Uncertain_significance'
  allowedClinVar['Uncertain_significance,_association']          = 'Uncertain_significance'
  allowedClinVar['Uncertain_significance|association']           = 'Uncertain_significance'
  allowedClinVar['Uncertain_significance|_association']          = 'Uncertain_significance'
  allowedClinVar['Uncertain_significance,_drug_response']        = 'Uncertain_significance'
  allowedClinVar['Uncertain_significance|drug_response']         = 'Uncertain_significance'
  allowedClinVar['Uncertain_significance|_drug_response']        = 'Uncertain_significance'
  allowedClinVar['Uncertain_significance,_Affects']              = 'Uncertain_significance'
  allowedClinVar['Uncertain_significance|Affects']               = 'Uncertain_significance'
  allowedClinVar['Uncertain_significance|_Affects']              = 'Uncertain_significance'
  allowedClinVar['Uncertain_significance|risk_factor']           = 'Uncertain_significance'
  allowedClinVar['Uncertain_significance/Uncertain_risk_allele'] = 'Uncertain_significance'

  allowedClinVar['Conflicting_interpretations_of_pathogenicity,_other']                    = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity|other']                     = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity|_other']                    = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity,_risk_factor']              = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity|risk_factor']               = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity|_risk_factor']              = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity,_protective']               = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity|protective']                = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity|_protective']               = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity,_Affects']                  = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity,_other,_risk_factor']       = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity|other|risk_factor']         = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity|_other|_risk_factor']       = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity,_drug_response']            = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity|drug_response']             = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity|_drug_response']            = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity,_association']              = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity|association']               = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity|_association']              = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity,_drug_response,_other']     = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity|drug_response|other']       = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity|_drug_response|_other']     = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity|association|risk_factor']   = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity,_association,_risk_factor'] = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity|_association|_risk_factor'] = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity|Affects']                   = 'Conflicting_interpretations_of_pathogenicity'
  allowedClinVar['Conflicting_interpretations_of_pathogenicity|_Affects']                  = 'Conflicting_interpretations_of_pathogenicity'

  allowedClinVar['Likely_benign,_other']                = 'Likely_benign'
  allowedClinVar['Likely_benign|other']                 = 'Likely_benign'
  allowedClinVar['Likely_benign|_other']                = 'Likely_benign'
  allowedClinVar['Likely_benign,_drug_response']        = 'Likely_benign'
  allowedClinVar['Likely_benign,_drug_response,_other'] = 'Likely_benign'
  allowedClinVar['Likely_benign|_drug_response|_other'] = 'Likely_benign'
  allowedClinVar['Likely_benign|drug_response|other']   = 'Likely_benign'
  allowedClinVar['Likely_benign|risk_factor']           = 'Likely_benign'
  allowedClinVar['Likely_benign|_risk_factor']          = 'Likely_benign'
  allowedClinVar['Likely_benign,_association']          = 'Likely_benign'
  allowedClinVar['Likely_benign|association' ]          = 'Likely_benign'

  allowedClinVar['Benign/Likely_benign,_risk_factor']          = 'Benign/Likely_benign'
  allowedClinVar['Benign/Likely_benign|risk_factor']           = 'Benign/Likely_benign'
  allowedClinVar['Benign/Likely_benign|_risk_factor']          = 'Benign/Likely_benign'
  allowedClinVar['Benign/Likely_benign,_other']                = 'Benign/Likely_benign'
  allowedClinVar['Benign/Likely_benign|other']                 = 'Benign/Likely_benign'
  allowedClinVar['Benign/Likely_benign|_other']                = 'Benign/Likely_benign'
  allowedClinVar['Benign/Likely_benign|other|risk_factor']     = 'Benign/Likely_benign'
  allowedClinVar['Benign/Likely_benign|_other|_risk_factor']   = 'Benign/Likely_benign'
  allowedClinVar['Benign/Likely_benign|association']           = 'Benign/Likely_benign'
  allowedClinVar['Benign/Likely_benign|_association']          = 'Benign/Likely_benign'
  allowedClinVar['Benign/Likely_benign,_drug_response']        = 'Benign/Likely_benign'
  allowedClinVar['Benign/Likely_benign|drug_response']         = 'Benign/Likely_benign'
  allowedClinVar['Benign/Likely_benign|_drug_response']        = 'Benign/Likely_benign'
  allowedClinVar['Benign/Likely_benign|drug_response|other']   = 'Benign/Likely_benign'
  allowedClinVar['Benign/Likely_benign|_drug_response|_other'] = 'Benign/Likely_benign'
  allowedClinVar['Benign/Likely_benign,_other,_risk_factor']   = 'Benign/Likely_benign'
  allowedClinVar['Benign/Likely_benign,_drug_response,_other'] = 'Benign/Likely_benign'

  allowedClinVar['Benign,_other']                            = 'Benign'
  allowedClinVar['Benign|other']                             = 'Benign'
  allowedClinVar['Benign|_other']                            = 'Benign'
  allowedClinVar['Benign,_risk_factor']                      = 'Benign'
  allowedClinVar['Benign|risk_factor']                       = 'Benign'
  allowedClinVar['Benign|_risk_factor']                      = 'Benign'
  allowedClinVar['Benign,_drug_response']                    = 'Benign'
  allowedClinVar['Benign|drug_response']                     = 'Benign'
  allowedClinVar['Benign|_drug_response']                    = 'Benign'
  allowedClinVar['Benign,_association']                      = 'Benign'
  allowedClinVar['Benign|association']                       = 'Benign'
  allowedClinVar['Benign|_association']                      = 'Benign'
  allowedClinVar['Benign,_protective']                       = 'Benign'
  allowedClinVar['Benign|protective']                        = 'Benign'
  allowedClinVar['Benign|_protective']                       = 'Benign'
  allowedClinVar['Benign|confers_sensitivity']               = 'Benign'
  allowedClinVar['Benign|_confers_sensitivity']              = 'Benign'
  allowedClinVar['Benign,_confers_sensitivity']              = 'Benign'
  allowedClinVar['Benign|_association|_confers_sensitivity'] = 'Benign'
  allowedClinVar['Benign,_association,_confers_sensitivity'] = 'Benign'
  allowedClinVar['Benign|association|confers_sensitivity']   = 'Benign'

  # Loop over all records in the vcf file
  for record in os.popen(bcftools.query(bcftoolsExe, vcf, tags)).readlines():
    noRecords += 1

    # Split the record on tabs
    fields = record.rstrip().split('\t')
  
    # Check that the variant has a value. If not, skip this record. The only records that will be output are those with values
    if fields[5] != '.':
  
      # Ensure the value is allowed. If the value is not recognized, store the value and skip the line
      if fields[5] not in allowedClinVar: 
        if fields[5] not in skippedValues: skippedValues[fields[5]] = 1
        else: skippedValues[fields[5]] += 1
        noSkippedRecords += 1
        continue
        #fail('Unknown ClinVar significance: "' + str(fields[5]) + '"')

      # Update the value if required
      else: fields[5] = allowedClinVar[fields[5]]

      # Update the chromosome and position
      fields[0], fields[2] = updateCoords(fields[0], fields[2])
  
      # Build the output record from the updated fields
      print('\t'.join(fields), file = outputFile)
      noOutputRecords += 1

    # Record the number of records with no value
    else: noDotRecords += 1

  # Write to screen all the records that were skipped based on unrecognised values
  if len(skippedValues.keys()) > 0:
    print('The following unrecognised values were observed resulting in the given number of skipped records:')
    for record in skippedValues: print('  ', record, ': ', skippedValues[record], sep = '')
  print()
  print('Total number of records in the input vcf file: ', noRecords, sep = '')
  print('Number of records with "." as the value      : ', noDotRecords, sep = '')
  print('Number of skipped records                    : ', noSkippedRecords, sep = '')
  print('Total number of output records               : ', noOutputRecords, sep = '')

# Update the chromosome and position in the tsv file
def updateCoords(chrom, pos):

  # Check that the chromosome does not include "chr" prefix
  if chrom.startswith('chr'): chrom = chrom.strip('chr')

  # Add one to the end position
  pos = str(int(pos) + 1)

  # Return the updated values
  return chrom, pos

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)

# Initialise global variables

# Pipeline version
version = "0.0.1"

# Define the allowed annotation classes
allowedClasses = ['A', 'B', 'C', 'clinvar', 'compound', 'OMIM', 'spliceai']

if __name__ == "__main__":
  main()
