import os

# Depending on the family type, certain aspects of the pipeline will change. Determine which sections of the
# pipeline can be omitted in this case
def determinePipeline(familyType):

  # Compound hets should only be included if both parents are present
  useInheritance = True if familyType == 'trio' else False

  # Define an object to contain all pipeline modifiers
  pipelineModifiers = {'useInheritance': useInheritance}

  # Return information on pipeline modifiers
  return pipelineModifiers

# Open a file to write the annotation pipeline script to
def openBashScript(workingDir):

  # Create a script file
  bashFilename  = workingDir + 'calypso_annotation_pipeline.sh'
  try: bashFile = open(bashFilename, "w")
  except: fail("There was a problem opening a file (calypso_annotation_pipeline.sh) to write to")

  # Return the file
  return bashFilename, bashFile

# Write all the resource file information to the bash file for the script to use
def bashResources(resourceInfo, workingDir, bashFile, vcf, chrFormat, ped, luaFilename, tomlFilename): #, familyType, pipelineModifiers):

  # Initial information
  #print('#! /bin/bash', file = bashFile)
  print('set -eou pipefail', file = bashFile)
  print(file = bashFile)

  # Define the tools to use
  print('# Define the tools to use while executing the Calypso pipeline', file = bashFile)
  print('DATAPATH=', resourceInfo['path'], sep = '', file = bashFile)
  if resourceInfo['toolsPath']:
    print('TOOLPATH=', resourceInfo['toolsPath'], sep = '', file = bashFile)
    print('BCFTOOLS=$TOOLPATH/bcftools/bcftools', sep = '', file = bashFile)
    print('export BCFTOOLS_PLUGINS=$TOOLPATH/bcftools/plugins', file = bashFile)
    print('VEP=$TOOLPATH/ensembl-vep/vep', sep = '', file = bashFile)
    print('SLIVAR=$TOOLPATH/slivar', sep = '', file = bashFile)
    print('VCFANNO=$TOOLPATH/vcfanno', sep = '', file = bashFile)
  else:
    print('BCFTOOLS=bcftools', file = bashFile)
    print('VEP=vep', file = bashFile)
    print('SLIVAR=slivar', file = bashFile)
    print('VCFANNO=vcfanno', file = bashFile)
  print(file = bashFile)

  # Define the names of the input and output files
  print('# Following are the input VCF and output files created by the pipeline', file = bashFile)
  print('VCF=', vcf, sep = '', file = bashFile)
  print(file = bashFile)

  # Generate the names of the intermediate and final vcf files
  vcfBase     = os.path.abspath(vcf).split('/')[-1].rstrip('vcf.gz')
  filteredVcf = str(vcfBase) + '_calypso_filtered.vcf.gz'
  slivarVcf   = str(vcfBase) + '_calypso_slivar.vcf.gz'
  comphetVcf  = str(vcfBase) + '_calypso_comphet.vcf'
  clinvarVcf  = str(vcfBase) + '_calypso_clinVar.vcf.gz'
  print('FILEPATH=', workingDir, sep = '', file = bashFile)
  print('ANNOTATEDVCF=$FILEPATH/' + str(vcfBase) + '_annotated.vcf.gz', sep = '', file = bashFile)

  # If compound hets are not generated, there is no COMPHETS output
  print('FILTEREDVCF=$FILEPATH/' + str(filteredVcf), sep = '', file = bashFile)
  print('SLIVARVCF=$FILEPATH/' + str(slivarVcf), sep = '', file = bashFile)
  print('COMPHETSVCF=$FILEPATH/' + str(comphetVcf), sep = '', file = bashFile)
  print('STDOUT=calypso_annotation_pipeline.stdout', file = bashFile)
  print('STDERR=calypso_annotation_pipeline.stderr', file = bashFile)

  # Write the ped file, if necessary
  print('PED=', ped, sep = '', file = bashFile)
  print(file = bashFile)

  # Define the required resources
  print('# Following is a list of required resources', file = bashFile)
  try: print('REF=$DATAPATH/', resourceInfo['resources']['fasta']['file'], sep = '', file = bashFile)
  except: fail('The resources json does not define a reference fasta file')
  try: print('GFF=$DATAPATH/', resourceInfo['resources']['gff']['file'], sep = '', file = bashFile)
  except: fail('The resources json does not define a gff file')
  print('TOML=$FILEPATH/', tomlFilename, sep = '', file = bashFile)
  print('LUA=$FILEPATH/', luaFilename, sep = '', file = bashFile)

  # If the chr map is required, include the path to it. The chromosome format is different depending on whether this is a GRCh37
  # or GRCh38 reference. GRCh37 needs to use the '1' format, but GRCh38 uses 'chr1'. chrFormat is False if in the '1' format
  if resourceInfo['reference'] == 'GRCh37':
    if chrFormat: 
      print('CHR_NOCHR_MAP=$DATAPATH/chr_nochr_map.txt', sep = '', file = bashFile)
      print('NOCHR_CHR_MAP=$DATAPATH/nochr_chr_map.txt', sep = '', file = bashFile)
  elif resourceInfo['reference'] == 'GRCh38':
    if not chrFormat:
      print('CHR_NOCHR_MAP=$DATAPATH/chr_nochr_map.txt', sep = '', file = bashFile)
      print('NOCHR_CHR_MAP=$DATAPATH/nochr_chr_map.txt', sep = '', file = bashFile)

  # The slivar js file is based on the family structure. Check that a file exists for the current family
  jsFile = 'slivar-functions.js'
  print('JS=$DATAPATH/slivar/', jsFile, sep = '', file = bashFile)
  print(file = bashFile)

  # Return the name of the filtered vcf file
  return filteredVcf

# Generate a text file containing all the samples
def samplesFile(bashFile):

  # Generate a samples text file from the ped file
  print('# Generate a text file containing all samples in the family', file = bashFile)
  print('echo -n "Creating text file of samples..." > $STDOUT', file = bashFile)
  print('echo -n "Creating text file of samples..."', file = bashFile)
  print('tail -n+2 $PED | cut -f 2 | sort -u > samples.txt', file = bashFile)
  print('echo "complete"', file = bashFile)
  print('echo "complete" >> $STDOUT', file = bashFile)
  print(file = bashFile)

# Annotate the vcf file using bcftools, vcfanno, and VEP
def annotateVcf(resourceInfo, bashFile, chrFormat, samples):

  # Generate a comma separated list of samples to extract from the vcf file
  sampleList = ''
  for sample in samples:
    if samples[sample]['vcf_sample_name']: sampleList += samples[sample]['vcf_sample_name'] + ','
  sampleList = sampleList.rstrip(',')

  # Generate the list of annotations to be output from VEP
  annotationString = ''
  for annotation in resourceInfo['resources']['vep']['fields']: annotationString += annotation + ','
  annotationString = annotationString.rstrip(',')

  # Print out status messages
  print('# Normalize, subset, and annotate original VCF', file = bashFile)
  print('echo -n "Subsetting, normalizing, and annotating input VCF..."', file = bashFile)

  # Include additional resources for VEP
  if 'export' in resourceInfo['resources']['vep']:
    for newPath in resourceInfo['resources']['vep']['export']: print(newPath, file = bashFile)

################
################
################ Limiting to autosome and X was included due to errors, but Y and MT need to be included
################
################
  # Start the script by extracting the required samples from the vcf, and limiting to the autosome and X chromosomes.
  # Ensure the correct format is used for the chromosomes. chrFormat is False if in the '1' format
  if resourceInfo['reference'] == 'GRCh37':
    if not chrFormat:
      chroms = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X'
      print('$BCFTOOLS view -a -c 1 -s "', sampleList, '" -r "', chroms, '" $VCF 2>> $STDERR \\', sep = '', file = bashFile)
      print('  | $BCFTOOLS annotate -x INFO --rename-chrs $NOCHR_CHR_MAP 2>> $STDERR \\', file = bashFile)
    else:
      chroms = 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX'
      print('$BCFTOOLS view -a -c 1 -s "', sampleList, '" -r "', chroms, '" $VCF 2>> $STDERR \\', sep = '', file = bashFile)
      print('  | $BCFTOOLS annotate -x INFO --rename-chrs $CHR_NOCHR_MAP 2>> $STDERR \\', file = bashFile)
  elif resourceInfo['reference'] == 'GRCh38':
    if not chrFormat:
      chroms = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X'
      print('$BCFTOOLS view -a -c 1 -s "', sampleList, '" -r "', chroms, '" $VCF 2>> $STDERR \\', sep = '', file = bashFile)
      print('  | $BCFTOOLS annotate -x INFO --rename-chrs $NOCHR_CHR_MAP 2>> $STDERR \\', file = bashFile)
    else:
      chroms = 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX'
      print('$BCFTOOLS view -a -c 1 -s "', sampleList, '" -r "', chroms, '" $VCF 2>> $STDERR \\', sep = '', file = bashFile)
      print('  | $BCFTOOLS annotate -x INFO 2>> $STDERR \\', file = bashFile)

  # Normalize and decompose the vcf, then extract the variants present in the required samples
  print('  | $BCFTOOLS norm -m - -w 10000 -f $REF 2> $STDERR \\', file = bashFile)
  print('  | $BCFTOOLS view -a -c 1 -s "', sampleList, '" 2>> $STDERR \\', sep = '', file = bashFile)

  # Annotate the vcf file using both bcftools csq and vcfanno
  print('  | $BCFTOOLS csq -f $REF --ncsq 40 -l -g $GFF 2>> $STDERR \\', file = bashFile)
  print('  | $VCFANNO -lua $LUA -p 16 $TOML /dev/stdin 2>> $STDERR \\', file = bashFile)

  # Annotate with VEP unless it is to be ignored
  if not resourceInfo['resources']['vep']['ignore']:
    print('  | $VEP \\', file = bashFile)
    print('    --assembly ', resourceInfo['reference'], ' \\', sep = '', file = bashFile)
    print('    --fasta $REF \\', file = bashFile)
    print('    --cache \\', file = bashFile)
    print('    --dir_cache ', resourceInfo['resources']['vep']['cache'], ' \\', sep = '', file = bashFile)
    print('    --dir_plugins ', resourceInfo['resources']['vep']['plugins'], ' \\', sep = '', file = bashFile)
    print('    --offline \\', file = bashFile)
    print('    --vcf \\', file = bashFile)
    print('    --force \\', file = bashFile)
    print('    --check_existing \\', file = bashFile)
    print('    --quiet \\', file = bashFile)
    print('    --fork 40 \\', file = bashFile)
    print('    --format vcf \\', file = bashFile)
    print('    --force_overwrite \\', file = bashFile)
  
    # Add additional options defined in the resource json
    if 'options' in resourceInfo['resources']['vep']:
      for vepOption in resourceInfo['resources']['vep']['options']: print('    ', vepOption, ' \\', sep = '', file = bashFile)

    # Include all plugins
    if 'active_plugins' in resourceInfo['resources']['vep']:
      for plugin in resourceInfo['resources']['vep']['active_plugins']:
        if 'files' in resourceInfo['resources']['vep']['active_plugins'][plugin]: print('    --plugin ', plugin, ',', resourceInfo['resources']['vep']['active_plugins'][plugin]['files'], ' \\', sep = '', file = bashFile)
        else: print('    --plugin ', plugin, ' \\', sep = '', file = bashFile)
    print('    --fields "', annotationString, '" \\', sep = '', file = bashFile)
    #print('    --fields "', resourceInfo['resources']['vep']['fields'], '" \\', sep = '', file = bashFile)
    #print('    --fields "SYMBOL,Feature,IMPACT,Consequence,HGVSc,HGVSp,MaxEntScan_diff,MaxEntScan_alt,GeneSplicer,five_prime_UTR_variant_annotation,SIFT,PolyPhen,PUBMED,LoFtool" \\', file = bashFile)
    print('    --no_stats \\', file = bashFile)
    print('    --output_file STDOUT \\', file = bashFile)
    print('    2>> $STDERR \\', file = bashFile)

  # Add VAF to the FORMAT and HWE to the info field
  print('  | $BCFTOOLS +fill-tags - -- -t FORMAT/VAF,HWE \\', file = bashFile)

  # Output the final vcf file. If the original file had the '1' format (not 'chr1') for chromosomes,
  # ensure the final vcd file is in this format
  if resourceInfo['reference'] == 'GRCh37':
    if chrFormat: print('  | $BCFTOOLS annotate -O z -o $ANNOTATEDVCF --rename-chrs $NOCHR_CHR_MAP - \\', file = bashFile)
    else: print('  | $BCFTOOLS view -O z -o $ANNOTATEDVCF - \\', file = bashFile)
  elif resourceInfo['reference'] == 'GRCh38':
    if not chrFormat: print('  | $BCFTOOLS annotate -O z -o $ANNOTATEDVCF --rename-chrs $CHR_NOCHR_MAP - \\', file = bashFile)
    else: print('  | $BCFTOOLS view -O z -o $ANNOTATEDVCF - \\', file = bashFile)
  else: fail('Unknown reference: ' + str(resourceInfo['reference']))
  print('  >> $STDOUT 2>> $STDERR', file = bashFile)
  print('echo "complete"', file = bashFile)
  print(file = bashFile)

# Filter the final vcf file to only include variants present in the proband, and based on some
# basic annotations
def filterVariants(bashFile, samples, proband, resourceInfo):

  # Create a file containing the name of the proband for use in the filter. Note that there can be multiple probands, e.g.
  # if the family is two affected siblings
  print('# Create file containing the proband(s) only', file = bashFile)
  print('echo -n "Creating proband(s) file..."', file = bashFile)
  print('touch proband.txt', file = bashFile)
  #for sample in proband: print('echo ', sample, ' >> proband.txt', sep = '', file = bashFile)
  print('echo ', samples[proband]['vcf_sample_name'], ' >> proband.txt', sep = '', file = bashFile)
#  for sample in samples:
#    print('TEST', sample, samples[sample])
#    if samples[sample]['isAffected']: print('echo ', sample, ' >> proband.txt', sep = '', file = bashFile)
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
  print('$BCFTOOLS view \\', file = bashFile)
  print('  -i \'GT[@proband.txt]="alt" & GT[@proband.txt]!="mis"\' \\', file = bashFile)
  print('  $ANNOTATEDVCF \\', file = bashFile)
  print('  2>> $STDERR \\', file = bashFile)
  print('  | $SLIVAR expr --vcf - \\', file = bashFile)

  # The filter should return all variants that are not in, or have popmax AF < 0.01 in gnomAD (v2.1.1 and v3.2.1)...
  print('  --info \'( ( (!("ge4_0_0_AF_grpmax" in INFO) || ("ge4_0_0_AF_grpmax" in INFO && INFO.ge4_0_0_AF_grpmax < 0.01)) && ', file = bashFile)
  print('            (  !("gg4_0_0_AF_grpmax" in INFO) || ("gg4_0_0_AF_grpmax" in INFO && INFO.gg4_0_0_AF_grpmax < 0.01)) ) ||', file = bashFile)

  # ...or that the variant is present in ClinVar...
  print('  ("CLNSIG" in INFO) ||', file = bashFile)

  # ...or that the CADD_Phred score is over 15...
  print('  ("CADD_Phred" in INFO && INFO.CADD_Phred > 15) ||', file = bashFile)

  # ...or that the REVEL or MutScores are over 0.1...
  print('  ("REVEL" in INFO && INFO.REVEL > 0.1) ||', file = bashFile)
  print('  ("MutScore" in INFO && INFO.MutScore > 0.1) ||', file = bashFile)

  # ...or that the variant has a spliceAI score > 0.1 for any of the acceptor / donor sites.
  print('  ("SpliceAI_DS_AG" in INFO && INFO.SpliceAI_DS_AG > 0.1) ||', file = bashFile)
  print('  ("SpliceAI_DS_AL" in INFO && INFO.SpliceAI_DS_AL > 0.1) ||', file = bashFile)
  print('  ("SpliceAI_DS_DG" in INFO && INFO.SpliceAI_DS_DG > 0.1) ||', file = bashFile)
  print('  ("SpliceAI_DS_DL" in INFO && INFO.SpliceAI_DS_DL > 0.1) ) &&', file = bashFile)

  # ...and that the variant does not have "*" as the alt allele
  print('  variant.FILTER == "PASS" && variant.ALT[0] != "*"\' \\', file = bashFile)
  print('  --pass-only \\', file = bashFile)
  print('  2>> $STDERR \\', file = bashFile)
#  print('  | $VEP \\', file = bashFile)
#  print('    --fasta $REF \\', file = bashFile)
#  print('    --dir_cache ', resourceInfo['resources']['vep']['cache'], ' \\', sep = '', file = bashFile)
#  print('    --dir_plugins ', resourceInfo['resources']['vep']['plugins'], ' \\', sep = '', file = bashFile)
#  print('    --assembly ', resourceInfo['reference'], '\\', sep = '', file = bashFile)
#  print('    --vcf \\', file = bashFile)
#  print('    --force \\', file = bashFile)
#  print('    --check_existing \\', file = bashFile)
#  print('    --quiet \\', file = bashFile)
#  print('    --fork 40 \\', file = bashFile)
#  print('    --format vcf \\', file = bashFile)
#  print('    --force_overwrite \\', file = bashFile)
#  print('    --cache \\', file = bashFile)
#  print('    --offline \\', file = bashFile)
#  print('    --no_stats \\', file = bashFile)
#  print('    --hgvs \\', file = bashFile)
#  print('    --plugin UTRannotator \\', file = bashFile)
#  print('    --fields "SYMBOL,Feature,IMPACT,Consequence,HGVSc,HGVSp" \\', file = bashFile)
#  print('    --output_file STDOUT \\', file = bashFile)
#  print('    2>> $STDERR \\', file = bashFile)
  print('  | $BCFTOOLS view -O z -o $FILTEREDVCF - \\', file = bashFile)
  print('  >> $STDOUT 2>> $STDERR', file = bashFile)
  print('echo "complete"', file = bashFile)
  print('$BCFTOOLS index -t $FILTEREDVCF', file = bashFile)
  print(file = bashFile)

# Use Slivar to extract compound heterozygous variants based on the Slivar rare disease wiki
def compHets(bashFile, rootPath, configFile, toolsDir, utilsDir, projectId):
  print('# Extract compound hets', file = bashFile)
  print('echo -n "Determining compound het variants..."', file = bashFile)
  print('$SLIVAR expr --vcf $FILTEREDVCF \\', file = bashFile)
  print('  --ped $PED \\', file = bashFile)
  print('  --pass-only \\', file = bashFile)
  print('  --js $JS \\', file = bashFile)
  print('  --info \'INFO.impactful && variant.FILTER == "PASS" && variant.ALT[0] != "*"\' \\', file = bashFile)
  print('  --family-expr \'denovo:fam.every(segregating_denovo)\' \\', file = bashFile)
  print('  --family-expr \'recessive:fam.every(segregating_recessive)\' \\', file = bashFile)
  print('  --family-expr \'x_denovo:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_denovo_x)\' \\', file = bashFile)
  print('  --family-expr \'x_recessive:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_recessive_x)\' \\', file = bashFile)
  print('  --trio \'comphet_side:comphet_side(kid, mom, dad)\' \\', file = bashFile)
  print('  -o $SLIVARVCF \\', file = bashFile)
  print('  >> $STDOUT 2>> $STDERR', file = bashFile)
  print('$BCFTOOLS index -t $SLIVARVCF', file = bashFile)
  print(file = bashFile)
  print('$SLIVAR compound-hets -v $SLIVARVCF \\', file = bashFile)
  print('  --sample-field comphet_side \\', file = bashFile)
  print('  --sample-field denovo \\', file = bashFile)
  print('  -p $PED > $COMPHETSVCF', file = bashFile)
  print('echo "complete"', file = bashFile)

  # The comp hets vcf contains all comp het variants. Get the gene names for these variants
  print(file = bashFile)
  print('# Extract all the gene names from the comp hets file', file = bashFile)
  print('echo -n "Extract all the gene names from the comp hets file..."', file = bashFile)
  print('python3 ', str(rootPath), '/scripts/comp_het_genes.py \\', sep = '', file = bashFile)
  print('  -c ', configFile, ' \\', sep = '', file = bashFile)
  print('  -s ', toolsDir, ' \\', sep= '', file = bashFile)
  print('  -l ', utilsDir, ' \\', sep= '', file = bashFile)
  print('  -p ', projectId, ' \\', sep = '', file = bashFile)
  print('  -i $COMPHETSVCF', sep = '', file = bashFile)
  print('echo "complete"', file = bashFile)
  print(file = bashFile)

# Delete files no longer required
def deleteFiles(args, pipelineModifiers, bashFile):
  print('# Delete files no longer required', file = bashFile)
  print('echo -n "Deleting files..."', file = bashFile)
  if pipelineModifiers['useInheritance']:
    print('rm -f $COMPHETSVCF', file = bashFile)
    print('rm -f $COMPHETSVCF.tbi', file = bashFile)
  print('rm -f samples.txt', file = bashFile)
  print('rm -f proband.txt', file = bashFile)
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
