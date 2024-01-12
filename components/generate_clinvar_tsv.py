from datetime import date
from os.path import exists
from sys import path

import argparse
import os
import sys
import calypso_resources as res

# Add the path of the common functions and import them
path.append(os.path.dirname(__file__) + '/components')
import tools_bcftools as bcftools
import calypso_mosaic_resources as mosr

def main():
  global allowedClinVar
  global bcftoolsExe

  # Parse the command line
  args = parseCommandLine()

  # Define the executable bcftools command
  if args.tools_directory:
    if not args.tools_directory.endswith('/'): args.tools_directory = str(args.tools_directory) + '/'
    bcftoolsExe = args.tools_directory + 'bcftools/bcftools'
  else: bcftoolsExe = 'bcftools'

  # Read the mosaicJson file to get information on how to process different annotations
  mosaicInfo = mosr.readMosaicJson(args.mosaic_json, args.reference)

  # Open an output tsv file to write annotations to
  outputFile = open(args.output_tsv, 'w')

  # Loop over the annotations that are to be uploaded to Mosaic for this resource and get the annotation name, uid and
  # type
  annotations = {}
  uids        = []
  for annotation in mosaicInfo['resources']['clinVar']['annotations']:
    uid     = mosaicInfo['resources']['clinVar']['annotations'][annotation]['uid']
    tagType = mosaicInfo['resources']['clinVar']['annotations'][annotation]['type']
    annotations[annotation] = {'uid': uid, 'type': tagType}

  # Write the header line to the tsv file
  print('CHROM\tSTART\tEND\tREF\tALT\t', '\t'.join(str(x['uid']) for x in annotations.values()), sep = '', file = outputFile)

  # Define some variables
  skippedValues    = {}
  noRecords        = 0
  noOutputRecords  = 0
  noSkippedRecords = 0
  noDotRecords     = 0

  # Loop over all records in the vcf file
  for record in os.popen(bcftools.query(bcftoolsExe, args.input_vcf, ['CLNSIG'])).readlines():
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
  
  # Close the output tsv file
  outputFile.close()

# Input options
def parseCommandLine():
  global version
  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  #parser.add_argument('--config', '-c', required = True, metavar = "string", help = "The config file for Mosaic")
  parser.add_argument('--input_vcf', '-i', required = True, metavar = 'string', help = 'The input vcf file to annotate')
  parser.add_argument('--output_tsv', '-o', required = True, metavar = 'string', help = 'The output tsv file')
  #parser.add_argument('--resource', '-e', required = True, metavar = 'string', help = 'The name of the resource (used for the output tsv name')
  parser.add_argument('--reference', '-r', required = True, metavar = 'string', help = 'The genome reference file used')
  #parser.add_argument('--tags', '-g', required = True, metavar = 'string', help = 'A comma separated list of VCF INFO tags to be extracted')
  #parser.add_argument('--uids', '-d', required = True, metavar = 'string', help = 'A comma separated list of uids for the resource annotations')
  #parser.add_argument('--utils_directory', '-l', required = True, metavar = 'string', help = 'The path to the public-utils directory')
  parser.add_argument('--tools_directory', '-s', required = True, metavar = 'string', help = 'The path to the directory where the tools live')
  #parser.add_argument('--resource_json', '-j', required = True, metavar = 'string', help = 'The json file describing the annotation resources')
  parser.add_argument('--mosaic_json', '-m', required = True, metavar = 'string', help = 'The json file describing the Mosaic parameters')

  # Optional mosaic arguments
  #parser.add_argument('--token', '-t', required = False, metavar = "string", help = "The Mosaic authorization token")
  #parser.add_argument('--url', '-u', required = False, metavar = "string", help = "The base url for Mosaic")
  #parser.add_argument('--attributes_project', '-a', required = False, metavar = "integer", help = "The Mosaic project id that contains public attributes")

  # Version
  parser.add_argument('--version', '-v', action="version", version='Calypso annotation pipeline version: ' + str(version))

  return parser.parse_args()

# The majority of annotations can be processed in a standard way. The type is checked and numerical
# annotations are checked to ensure they fall within the required bounds
def processStandard(vcf, annotations, outputFile):
  global bcftoolsExe

  # Loop over all records in the vcf file
  for record in os.popen(bcftools.query(bcftoolsExe, vcf, annotations.keys())).readlines():

    # Split the record on tabs
    fields = record.rstrip().split('\t')

    # If all values are '.', this line can be ignored
    uniqueValues = set(fields[5:])
    if len(uniqueValues) == 1 and list(uniqueValues)[0] == '.': continue

    # Update the chromosome and position
    fields[0], fields[2] = updateCoords(fields[0], fields[2])
  
    # Loop over the annotations in the order they appear in the bcftools query command
    i = 5
    for annotation in annotations:

      # If the annotation is numeric
      if str(annotations[annotation]['type']) == 'integer' or str(annotations[annotation]['type']) == 'float':

        # If the are multiple values (e.g. there is a ',') in the value, ensure that all values are within
        # the limits (1E-37 < x < 1E37), but output all values
        values      = fields[i].split(',') if ',' in fields[i] else [fields[i]]
        finalValues = []
        for value in values:
          if value == '.': finalValues.append('')
          else:
            try: typeTest = float(value)
            except: fail('Invalid value for annotation: ' + record.rstrip())
            if abs(float(value)) <= 1e-37: finalValues.append('0')
            elif abs(float(value)) >= 1e37: finalValues.append('1e37')
            else: finalValues.append(value)

        # Join the values back together with commas
        fields[i] = ','.join(finalValues)

      # Else if the annotation is a string, the annotation should be output as is
      elif str(annotations[annotation]['type']) == 'string': pass

      # If the type is unknown, fail
      else: fail('Annotation type (' + str(annotations[annotation]['type']) + ') for annotation ' + str(annotation) + ', is unknown')

      # Iterate to the next annotation field
      i += 1

    # Build the output record from the updated fields
    print('\t'.join(fields), file = outputFile)
  
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

# Define the bcftools executable
bcftoolsExe = False

# Define the allowed values for the CLNSIG field
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
allowedClinVar['Pathogenic/Likely_risk_allele|risk_factor']          = 'Pathogenic'
allowedClinVar['Pathogenic/Pathogenic,_low_penetrance|other']        = 'Pathogenic'

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
allowedClinVar['Likely_pathogenic|protective']          = 'Likely_pathogenic'

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
allowedClinVar['Benign|Affects']                           = 'Benign'

if __name__ == "__main__":
  main()
