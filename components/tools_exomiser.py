#!/usr/bin/python

from os.path import exists

import os

########
########
######## This should be updated to get versions from json resource file
########
########
# Generate the application properties file
def applicationProperties(workingDir, toolsDir):

  # Create the application.properties file
  propFileName = str(workingDir) + 'application.properties'
  propFile     = open(propFileName, 'w')

  # Write the data versions to the file
  print('exomiser.data-directory=', toolsDir, 'exomiser-cli-13.3.0/data', sep = '', file = propFile)
  print('remm.version=0.3.1.post1', file = propFile)
  print('cadd.version=1.4', file = propFile)
  print('exomiser.hg38.data-version=2302', file = propFile)
  print('exomiser.phenotype.data-version=2302', file = propFile)

  # Close the application properties file
  propFile.close()

# Generate the yml file for running exomiser
def generateYml(workingDir, proband, reference, vcf, ped, hpo):

  # Determine if exomiser should use resources for hg19 or hg38
  if reference == 'GRCh37': exomiserRef = 'hg19'
  elif reference == 'GRCh38': exomiserRef = 'hg38'
  else: fail('Unknown reference for use with exomiser: ' + str(reference))

  # Generate a list of hpo terms, with the terms contained in single quotes
  hpoString = '\', \''.join(hpo)

  # Open a new yaml file
  yamlName = 'exomiser_' + str(proband) + '.yml'
  yaml     = open(str(workingDir) + yamlName, 'w')

  # Write the information about this project to the yaml
  print('# Exomiser analysis file for proband: ', proband, sep = '', file = yaml)
  print(file = yaml)
  print('analysis:', file = yaml)
  print('  genomeAssembly: ', str(exomiserRef), sep = '', file = yaml)
  print('  vcf: ', vcf, sep = '', file = yaml)
  print('  ped: ', ped, sep = '', file = yaml)
  print('  proband: ', proband, sep = '', file = yaml)
  print('  hpoIds: [\'', hpoString, '\']', sep = '', file = yaml)
  print(file = yaml)

  # Include the inheritance modes to consider
  # These are the default settings, with values representing the maximum minor allele frequency in percent (%) permitted for an
  # allele to be considered as a causative candidate under that mode of inheritance.
  # If you just want to analyse a sample under a single inheritance mode, delete/comment-out the others. For AUTOSOMAL_RECESSIVE
  # or X_RECESSIVE ensure *both* relevant HOM_ALT and COMP_HET modes are present.
  # In cases where you do not want any cut-offs applied an empty map should be used e.g. inheritanceModes: {}
  print('  inheritanceModes: {', file = yaml)
  print('    AUTOSOMAL_DOMINANT: 0.1,', file = yaml)
  print('    AUTOSOMAL_RECESSIVE_HOM_ALT: 0.1,', file = yaml)
  print('    AUTOSOMAL_RECESSIVE_COMP_HET: 2.0,', file = yaml)
  print('    X_DOMINANT: 0.1,', file = yaml)
  print('    X_RECESSIVE_HOM_ALT: 0.1,', file = yaml)
  print('    X_RECESSIVE_COMP_HET: 2.0,', file = yaml)
  print('    MITOCHONDRIAL: 0.2', file = yaml)
  print('  }', file = yaml)
  print(file = yaml)

  # Use the PASS_ONLY analysis mode
  print('  analysisMode: PASS_ONLY', file = yaml)
  print(file = yaml)

  # Define the allele frequency sources to be used
  print('  frequencySources: [', file = yaml)
  print('    THOUSAND_GENOMES,', file = yaml)
  print('    TOPMED,', file = yaml)
  print('    UK10K,', file = yaml)
  print('    ESP_AFRICAN_AMERICAN,', file = yaml)
  print('    ESP_EUROPEAN_AMERICAN,', file = yaml)
  print('    ESP_ALL,', file = yaml)
  print('    EXAC_AFRICAN_INC_AFRICAN_AMERICAN,', file = yaml)
  print('    EXAC_AMERICAN,', file = yaml)
  print('    EXAC_SOUTH_ASIAN,', file = yaml)
  print('    EXAC_EAST_ASIAN,', file = yaml)
  print('    EXAC_FINNISH,', file = yaml)
  print('    EXAC_NON_FINNISH_EUROPEAN,', file = yaml)
  print('    EXAC_OTHER,', file = yaml)
  print('    GNOMAD_E_AFR,', file = yaml)
  print('    GNOMAD_E_AMR,', file = yaml)
  print('    GNOMAD_E_EAS,', file = yaml)
  print('    GNOMAD_E_FIN,', file = yaml)
  print('    GNOMAD_E_NFE,', file = yaml)
  print('    GNOMAD_E_SAS,', file = yaml)
  print('    GNOMAD_E_OTH,', file = yaml)
  print('    GNOMAD_G_AFR,', file = yaml)
  print('    GNOMAD_G_AMR,', file = yaml)
  print('    GNOMAD_G_EAS,', file = yaml)
  print('    GNOMAD_G_FIN,', file = yaml)
  print('    GNOMAD_G_NFE,', file = yaml)
  print('    GNOMAD_G_SAS,', file = yaml)
  print('    GNOMAD_G_OTH', file = yaml)
  print('  ]', file = yaml)
  print(file = yaml)

  # Define the resources to be used for determining pathogenicity
  print('  pathogenicitySources: [', file = yaml)
  print('    REVEL,', file = yaml)
  print('    MVP', file = yaml)
  print('  ]', file = yaml)
  print(file = yaml)

  # Define the order of steps for exomiser to take
  print('  steps: [', file = yaml)
  print('    failedVariantFilter: { },', file = yaml)
  print('    variantEffectFilter: {', file = yaml)
  print('      remove: [', file = yaml)
  print('        FIVE_PRIME_UTR_EXON_VARIANT,', file = yaml)
  print('        FIVE_PRIME_UTR_INTRON_VARIANT,', file = yaml)
  print('        THREE_PRIME_UTR_EXON_VARIANT,', file = yaml)
  print('        THREE_PRIME_UTR_INTRON_VARIANT,', file = yaml)
  print('        NON_CODING_TRANSCRIPT_EXON_VARIANT,', file = yaml)
  print('        NON_CODING_TRANSCRIPT_INTRON_VARIANT,', file = yaml)
  print('        CODING_TRANSCRIPT_INTRON_VARIANT,', file = yaml)
  print('        UPSTREAM_GENE_VARIANT,', file = yaml)
  print('        DOWNSTREAM_GENE_VARIANT,', file = yaml)
  print('        INTERGENIC_VARIANT,', file = yaml)
  print('        REGULATORY_REGION_VARIANT', file = yaml)
  print('      ]', file = yaml)
  print('    },', file = yaml)
  print('    frequencyFilter: {maxFrequency: 2.0},', file = yaml)
  print('    pathogenicityFilter: {keepNonPathogenic: true},', file = yaml)
  print('    inheritanceFilter: {},', file = yaml)
  print('    omimPrioritiser: {},', file = yaml)
  print('    hiPhivePrioritiser: {runParams: \'human\'},', file = yaml)
  print('  ]', file = yaml)
  print(file = yaml)

  # Define the output options
  print('outputOptions:', file = yaml)
  print('  outputContributingVariantsOnly: false', file = yaml)
  print('  numGenes: 0', file = yaml)
  print('  outputDirectory: ', str(workingDir) + 'exomiser_results', file = yaml)
  print('  outputFileName: ', str(proband), file = yaml)
  print('  outputFormats: [HTML, JSON, TSV_GENE, TSV_VARIANT, VCF]', file = yaml)

  # Close the yaml file
  yaml.close()

  # Return the path and name of the yaml file
  return yamlName

# Generate a script file to run exomiser
def generateScript(workingDir, toolsDir, yaml):

  # Create script file for running exomiser
  scriptName = str(workingDir) + 'calypso_exomiser.sh'
  script     = open(scriptName, 'w')
  print('set -eou pipefail', file = script)
  print(file = script)

  # Create variables for paths
  print('TOOLSPATH=', toolsDir, sep = '', file = script)
  print('WORKINGPATH=', workingDir, sep = '', file = script)
  print('STDOUT=', workingDir, 'exomiser.stdout', sep = '', file = script)
  print('STDERR=', workingDir, 'exomiser.stderr', sep = '', file = script)
  print(file = script)

#####
##### UTAH SPECIFIC REQUIREMENT
#####
  print('# Module requirement for Utah environment', file = script)
  print('module load jdk/11', file = script)
  print(file = script)

#########
#########
######### Update with version from resources json
#########
#########
  # Include the executable for running exomiser
  print('java -jar $TOOLSPATH/exomiser-cli-13.3.0/exomiser-cli-13.3.0.jar \\', sep = '', file = script)
  print('  --analysis $WORKINGPATH/', yaml, ' \\', sep = '', file = script)
  print(' > $STDOUT \\', file = script)
  print(' 2> $STDERR', file = script)

  # Return the script
  return scriptName, script

# Parse the output variants.tsv file for upload to Mosaic
def parseOutput(script, scriptDir, config, utilsDir, probandName, projectId):

  # Define the name of the output tsv
  outputName = 'exomiser.tsv'

  # The pvalue for defining the top candidates is hard-coded to 0.01
  pvalue = 0.01

  print(file = script)
  print('# Parse the resulting output file', file = script)
  print('UTILSPATH=', str(utilsDir), sep = '', file = script)
  print(file = script)
  print('echo -n "Creating tsv file for exomiser variants..."', file = script)
  print('python ', str(scriptDir), '/parse_exomiser_variants.py \\', sep = '', file = script)
  print('  -c ', str(config), ' \\', sep = '', file = script)
  print('  -l ', str(utilsDir), ' \\', sep = '', file = script)
  print('  -i $WORKINGPATH/exomiser_results/', str(probandName), '.variants.tsv \\', sep = '', file = script)
  print('  -e ', str(pvalue), ' \\', sep = '', file = script)
  print('  -o  ', str(outputName), ' \\', sep = '', file = script)
  print('  -p ', str(projectId), sep = '', file = script)
  print('echo "complete"', file = script)

  # If the parsing completed successfully, upload the annotations to mosaic
  print(file = script)
  print('if [[ $? == 0 ]];', file = script)
  print('then', file = script)
  print('  echo -n "Uploading exomiser annotations..."', file = script)
  print('  python $UTILSPATH/scripts/upload_annotations.py \\', file = script)
  print('    -c ', str(config), ' \\', sep = '', file = script)
  print('    -i ', str(outputName), ' \\', sep = '', file = script)
  print('    -p ', str(projectId), sep = '', file = script)
  print('  echo "complete"', file = script)
  print('fi', file = script)

# Close the exomiser script and make it executable
def closeFile(scriptName, script):

  # Close the exomiser script file
  script.close()

  # Make the annotation script executable
  makeExecutable = os.popen('chmod +x ' + scriptName).read()

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
