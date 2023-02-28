#!/usr/bin/python

from __future__ import print_function

import os

# Open a file to write the annotation pipeline script to
def openBashScript(workingDir):

  # Create a script file
  bashFilename  = workingDir + 'calypso_annotation_pipeline.sh'
  try: bashFile = open(bashFilename, "w")
  except: fail("There was a problem opening a file (calypso_annotation_pipeline.sh) to write to")

  # Return the file
  return bashFilename, bashFile

# Write all the resource file information to the bash file for the script to use
def bashResources(resourceInfo, workingDir, bashFile, vcf, ped, tomlFilename):

  # Initial information
  print('#! /bin/bash', file = bashFile)
  print('set -eou pipefail', file = bashFile)
  print(file = bashFile)

  # Define the names of the input and output files
  print('# Following are the input VCF and output files created by the pipeline', file = bashFile)
  print('VCF=', os.path.abspath(vcf), sep = '', file = bashFile)

  # Generate the names of the intermediate and final vcf files
  vcfBase     = workingDir + os.path.abspath(vcf).split('/')[-1].rstrip('vcf.gz')
  filteredVcf = str(vcfBase) + '_calypso_filtered.vcf.gz'
  clinvarVcf  = str(vcfBase) + '_calypso_clinVar.vcf.gz'
  rareVcf     = str(vcfBase) + '_calypso_rare_disease.vcf.gz'
  print('CLEANVCF=' + str(vcfBase) + '_clean.vcf.gz', sep = '', file = bashFile)
  print('ANNOTATEDVCF=' + str(vcfBase) + '_annotated.vcf.gz', sep = '', file = bashFile)
  print('PROBANDVCF=' + str(vcfBase) + '_proband.vcf.gz', sep = '', file = bashFile)
  print('PREANNOVCF=' + str(vcfBase) + '_preanno.vcf.gz', sep = '', file = bashFile)
  print('COMPHETS=' + str(vcfBase) + '_comphets.vcf.gz', sep = '', file = bashFile)
  print('FINALVCF=' + str(vcfBase) + '_calypso.vcf.gz', sep = '', file = bashFile)
  print('FILTEREDVCF=' + str(filteredVcf), sep = '', file = bashFile)
  print('CLINVARVCF=' + str(clinvarVcf), sep = '', file = bashFile)
  print('RAREDISEASEVCF=' + str(rareVcf), sep = '', file = bashFile)
  print('STDOUT=calypso_annotation_pipeline.stdout', file = bashFile)
  print('STDERR=calypso_annotation_pipeline.stderr', file = bashFile)

  # Write the ped file, if necessary
  print('PED=', os.path.abspath(ped), sep = '', file = bashFile)
  print(file = bashFile)

  # Define the required resources
  print('# Following is a list of required resources', file = bashFile)
  try: print('REF=', resourceInfo['resources']['fasta']['file'], sep = '', file = bashFile)
  except: fail('The resources json does not define a reference fasta file')
  try: print('GFF=', resourceInfo['resources']['gff']['file'], sep = '', file = bashFile)
  except: fail('The resources json does not define a gff file')
  try: print('SLIVAR_GNOMAD=', resourceInfo['resources']['slivar_gnomAD']['file'], sep = '', file = bashFile)
  except: fail('The resources json does not define a gnomAD zip file')
  try: print('JS=', resourceInfo['resources']['slivar_js']['file'], sep = '', file = bashFile)
  except: fail('The resources json does not define the Slivar functions js file')
  print('TOML=', tomlFilename, sep = '', file = bashFile)
  print(file = bashFile)

  # Return the name of the filtered vcf file
  return filteredVcf, clinvarVcf, rareVcf

# Generate a text file containing all the samples
def samplesFile(bashFile):
  print('# Generate a text file containing all samples in the family', file = bashFile)
  print('echo -n "Creating text file of samples..." > $STDOUT', file = bashFile)
  print('echo -n "Creating text file of samples..."', file = bashFile)
  print('tail -n+2 $PED | cut -f 2 | sort -u > samples.txt', file = bashFile)
  print('echo "complete"', file = bashFile)
  print('echo "complete" >> $STDOUT', file = bashFile)
  print(file = bashFile)

# Generate a normalized vcf containing only family members and no annotations
def cleanedVcf(bashFile):

  # Print out status messages
  print('# Normalize and subset original VCF', file = bashFile)
  print('echo -n "Subsetting and normalizing input VCF..."', file = bashFile)

  # Generate a samples text file from the ped file
  print('bcftools norm -m - -w 10000 -f $REF $VCF 2> $STDERR\\', file = bashFile)
  print('  | bcftools view -a -c 1 -S samples.txt 2>> $STDERR \\', file = bashFile)
  print('  | bcftools annotate -x INFO -O z -o $CLEANVCF \\', file = bashFile)
  print('  > $STDOUT 2>> $STDERR', file = bashFile)
  print('echo "complete"', file = bashFile)
  print(file = bashFile)

  # Index the clean vcf file
  print('# Index the clean VCF file', file = bashFile)
  print('echo -n "Indexing cleaned VCF..."', file = bashFile)
  print('bcftools index -t $CLEANVCF', file = bashFile)
  print('echo "complete"', file = bashFile)
  print(file = bashFile)

# Annotate the cleaned vcf file using vcfanno
def annotateVcf(bashFile, resourceInfo, familyType):

  # Annotate the vcf with vcfanno
  print('# Annotate with bcftools csq and add further annotations with vcfanno', file = bashFile)
  print('echo -n "Annotating cleaned VCF..."', file = bashFile)
  print('bcftools csq -f $REF --ncsq 40 -l -g $GFF $CLEANVCF 2>> $STDERR \\', file = bashFile)
  print('  | vcfanno -p 16 $TOML /dev/stdin 2>> $STDERR \\', file = bashFile)

  # If this is a family, also include modes of inheritance using Slivar
  if familyType != 'singleton':
    print('  | slivar_static expr \\', file = bashFile)
    print('  --vcf /dev/stdin \\', file = bashFile)
    print('  --ped $PED \\', file = bashFile)
    print('  --js $JS \\', file = bashFile)
    print('  -g $SLIVAR_GNOMAD \\', file = bashFile)
    print('  --info \'INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"\' \\', file = bashFile)
    print('  --family-expr \'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < 0.001\' \\', file = bashFile)
    print('  --family-expr \'x_denovo:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af < 0.001\' \\', file = bashFile)
    print('  --family-expr \'recessive:fam.every(segregating_recessive)\' \\', file = bashFile)
    print('  --family-expr \'dominant:fam.every(segregating_dominant)\' \\', file = bashFile)
    print('  -o $ANNOTATEDVCF \\', file = bashFile)

  # If this is a singleton, the inheritance annoatations will not be performed and this is the final vcf
  else: print('  | bcftools view -O z -o $FINALVCF - \\', file = bashFile)
  print('  >> $STDOUT 2>> $STDERR', file = bashFile)
  print('echo "complete"', file = bashFile)
  print(file = bashFile)

  # Annotate with compound hets
  if familyType != 'singleton':

    # Annotate with compound hets
    print('# Annotate compound hets with Slivar', file = bashFile)
    print('echo -n "Adding in compound hets..."', file = bashFile)
    print('slivar_static expr \\', file = bashFile)
    print('  --pass-only \\', file = bashFile)
    print('  --vcf $ANNOTATEDVCF \\', file = bashFile)
    print('  --ped $PED \\', file = bashFile)
    print('  --js $JS \\', file = bashFile)
    print('  -g $SLIVAR_GNOMAD \\', file = bashFile)
    print('  --family-expr \'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < 0.001\' \\', file = bashFile)
    print('  --trio \'comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_popmax_af < 0.005\' \\', file = bashFile)
    print('  2>> $STDERR \\', file = bashFile)
    print('  | slivar_static compound-hets \\', file = bashFile)
    print('  --vcf /dev/stdin \\', file = bashFile)
    print('  --skip NONE \\', file = bashFile)
    print('  -s comphet_side \\', file = bashFile)
    print('  -s denovo \\', file = bashFile)
    print('  -p $PED \\', file = bashFile)
    print('  -o $COMPHETS \\', file = bashFile)
    print('  >> $STDOUT 2>> $STDERR', file = bashFile)
    print('echo "complete"', file = bashFile)
    print(file = bashFile)

    # Index the annotated VCF
    print('# Index the annotated VCF file', file = bashFile)
    print('echo -n "Indexing annotated VCF..."', file = bashFile)
    print('bcftools index -t $ANNOTATEDVCF', file = bashFile)
    print('echo "complete"', file = bashFile)
    print(file = bashFile)

    # Index the comphets vcf file
    print('# Index the comphets VCF file', file = bashFile)
    print('echo -n "Indexing comphets VCF..."', file = bashFile)
    print('bcftools index -t $COMPHETS', file = bashFile)
    print('echo "complete"', file = bashFile)
    print(file = bashFile)

    # Merge the compound hets into the annotated vcf
    print('# Merge annotated VCF with comphets', file = bashFile)
    print('echo -n "Merging annotated and compound hets VCF files..."', file = bashFile)
    print('bcftools annotate \\', file = bashFile)
    print('  -a $COMPHETS \\', file = bashFile)
    print('  -c \'INFO/slivar_comphet,INFO/comphet_side\' \\', file = bashFile)
    print('  -O z \\', file = bashFile)
    print('  -o $FINALVCF \\', file = bashFile)
    print('  $ANNOTATEDVCF \\', file = bashFile)
    print('  >> $STDOUT 2>> $STDERR', file = bashFile)
    print('bcftools index -t $FINALVCF', file = bashFile)
    print('echo "complete"', file = bashFile)
    print(file = bashFile)

# Filter the final vcf file to only include variants present in the proband, and based on some
# basic annotations
def filterVariants(bashFile, proband, resourceInfo):

  # Create a file containing the name of the proband for use in the filter
  print('# Create file containing the proband only', file = bashFile)
  print('echo -n "Creating proband file..."', file = bashFile)
  print('echo ', proband, ' > proband.txt', sep = '', file = bashFile)
  print('echo "complete"', file = bashFile)
  print(file = bashFile)

  # Filter the VCF file to generate variants to pass to Mosaic. This is a relatively permissive filter
  # to allow researchers to search variants. Annotate these variants with additional information with VEP
  #   1. Have PASS in the FILTER field
  #   2. Do not have * as the ALT allele
  #   3. Have a popmax AF < 0.01
  #   4. The proband has a genotype containing the alt allele
  print('# Filter the VCF file to generate variants to pass to Mosaic', file = bashFile)
  print('echo -n "Filtering final VCF file..."', file = bashFile)
  print('bcftools view \\', file = bashFile)
  print('  -i \'GT[@proband.txt]="alt" & GT[@proband.txt]!="mis"\' \\', file = bashFile)
  print('  $FINALVCF \\', file = bashFile)
  print('  2>> $STDERR \\', file = bashFile)
  print('  | slivar_static expr --vcf - \\', file = bashFile)
  print('  --info \'INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"\' \\', file = bashFile)
  print('  --pass-only \\', file = bashFile)
  print('  2>> $STDERR \\', file = bashFile)
  print('  | vep \\', file = bashFile)
  print('    --vcf \\', file = bashFile)
  print('    --force \\', file = bashFile)
  print('    --check_existing \\', file = bashFile)
  print('    --quiet \\', file = bashFile)
  print('    --fork 40 \\', file = bashFile)
  print('    --format vcf \\', file = bashFile)
  print('    --force_overwrite \\', file = bashFile)
  print('    --cache \\', file = bashFile)
  print('    --no_stats \\', file = bashFile)
  print('    --fasta $REF \\', file = bashFile)
  print('    --dir_cache ', resourceInfo['resources']['vep']['cache'], ' \\', sep = '', file = bashFile)
  print('    --dir_plugins ', resourceInfo['resources']['vep']['plugins'], ' \\', sep = '', file = bashFile)
  print('    --assembly ', resourceInfo['reference'], '\\', sep = '', file = bashFile)
  print('    --hgvs \\', file = bashFile)
  print('    --fields "SYMBOL,Feature,IMPACT,Consequence,HGVSc,HGVSp" \\', file = bashFile)
  print('    --output_file STDOUT \\', file = bashFile)
  print('    2>> $STDERR \\', file = bashFile)
  print('  | bcftools view -O z -o $FILTEREDVCF - \\', file = bashFile)
  print('  >> $STDOUT 2>> $STDERR', file = bashFile)
  print('echo "complete"', file = bashFile)
  print('bcftools index -t $FILTEREDVCF', file = bashFile)
  print(file = bashFile)

# Extract all variants that have "athogenic" in the ClinVar significance. This will extract all Pathogenic,
# Likely_pathogenic, and Conflicting... variants
def clinVarVariants(bashFile):
  print('# Extract all variants with some pathogenic significance', file = bashFile)
  print('echo -n "Generating clinVar variants file..."', file = bashFile)
  print('bcftools view -O z \\', file = bashFile)
  print('  -o $CLINVARVCF \\', file = bashFile)
  print('  -i \'INFO/CLNSIG~"athogenic"\' \\', file = bashFile)
  print('  $FINALVCF \\', file = bashFile)
  print('  >> $STDOUT 2>> $STDERR', file = bashFile)
  print('echo "complete"', file = bashFile)
  print('bcftools index -t $CLINVARVCF', file = bashFile)
  print(file = bashFile)

# Use Slivar to extract variants based on the Slivar rare disease wiki
def rareDiseaseVariants(bashFile):
  print('# Generate rare disease variants based on Slivar wiki', file = bashFile)
  print('echo -n "Generating rare disease variants..."', file = bashFile)
  print('bcftools csq -s - --ncsq 40 -g $GFF -l -f $REF $FILTEREDVCF -O u \\', file = bashFile)
  print('  2>> $STDERR \\', file = bashFile)
  print('  | slivar_static expr --vcf - \\', file = bashFile)
  print('  --ped $PED \\', file = bashFile)
  print('  -o $RAREDISEASEVCF \\', file = bashFile)
  print('  --pass-only \\', file = bashFile)
  print('  -g $SLIVAR_GNOMAD \\', file = bashFile)
  print('  --info \'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"\' \\', file = bashFile)
  print('  --js $JS \\', file = bashFile)
  print('  --family-expr \'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < 0.001\' \\', file = bashFile)
  print('  --family-expr \'recessive:fam.every(segregating_recessive)\' \\', file = bashFile)
  print('  --family-expr \'x_denovo:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af < 0.001\' \\', file = bashFile)
  print('  --family-expr \'x_recessive:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_recessive_x)\' \\', file = bashFile)
  print('  --trio \'comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_nhomalt < 10\' \\', file = bashFile)
  print('  >> $STDOUT 2>> $STDERR', file = bashFile)
  print('echo "complete"', file = bashFile)
  print('bcftools index -t $RAREDISEASEVCF', file = bashFile)
  print(file = bashFile)

# Delete files no longer required
def deleteFiles(args, bashFile, deletePed):
  print('# Delete files no longer required', file = bashFile)
  print('echo -n "Deleting files..."', file = bashFile)
  print('rm -f $CLEANVCF', file = bashFile)
  print('rm -f $CLEANVCF.tbi', file = bashFile)
  print('rm -f $ANNOTATEDVCF', file = bashFile)
  print('rm -f $ANNOTATEDVCF.tbi', file = bashFile)
  print('rm -f $COMPHETS', file = bashFile)
  print('rm -f $COMPHETS.tbi', file = bashFile)
  print('rm -f samples.txt', file = bashFile)
  print('rm -f proband.txt', file = bashFile)
  if deletePed: print('rm -f ', args.ped, sep = '', file = bashFile)
  print('echo "complete"', file = bashFile)
  print(file = bashFile)

# Close the script and make it executable
def finishScript(bashFile, bashFilename, version):

  # Print the pipeline completed successfully
  print('echo "Calypso pipeline version ', version, ' completed successfully"', sep = '', file = bashFile)

  # Close the file
  bashFile.close()

  # Make the annotation script executable
  makeExecutable = os.popen('chmod +x ' + bashFilename).read()

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
