import os
import json
import subprocess
import tools_bcftools as bcf

from os.path import exists

# If a vcf file was supplied on the command, get the sample name from the vcf header
def parseVcf(bcftools, vcf, mosaicSamples):

  # Check the file exists
  if not exists(vcf): fail('Input vcf file does not exist')

  # Get the vcf header
  header = os.popen(bcf.getHeader(bcftools, vcf)).read()
  for line in header.split('\n'):
    if line.startswith('#CHROM'):
      samples = line.rstrip().split('\t')[9:]
      break

  # Loop over the samples in the vcf file and find the ones in the mosaicSamples list
  for sample in samples:

#####
#####
##### UDN TWEAK TO BE FIXED. REMOVE the trailing "-XXXX" from the sample name in the comparison
#####
#####
    if '-' in sample: sample = sample.split('-')[0]
    if sample in mosaicSamples:
      mosaicSamples[sample]['vcf_file']        = vcf
      mosaicSamples[sample]['vcf_sample_name'] = sample
      print('  Sample ', sample, ' appears as "', sample, '" in the header of vcf file: ', vcf, sep = '')

  # Return the updated mosaicSamples and the vcf file
  return mosaicSamples, vcf

# Get the vcf files for all the mosaic samples
def getVcfFiles(mosaicConfig, api_sf, projectId, mosaicSamples):
  allVcfs = []

  # Loop over all the samples
  for sample in mosaicSamples:
    sampleId = mosaicSamples[sample]['id']
    vcfFiles = api_sf.getSampleVcfs(mosaicConfig, projectId, sampleId)

    # If there are no vcf files for this sample output a warning, but continue, unless this is the proband
    if len(vcfFiles) == 0:
      if mosaicSamples[sample]['relation'] == 'Proband': fail('calypso_vcf_files: The proband (' + sample + ') is not associated with any vcf files. The proband must appear in a vcf file (a vcf file can be specified on the command line using the --input_vcf (-i) argument)')
      else:
        print('  WARNING: Sample ' + sample + ', listed as having the relationship ' + mosaicSamples[sample]['relation'] + ' is not associated with a vcf and so will contribute no variants')
        mosaicSamples[sample]['vcf_file']        = False
        mosaicSamples[sample]['vcf_sample_name'] = False

    # If there is more than one vcf file for the sample, fail
    elif len(vcfFiles) > 1: fail('calypso_vcf_files: Mosaic has no, or multiple vcf files associated with it, so the file to use cannot be determined')
    else:
      for vcfFile in vcfFiles:
        uri = vcfFiles[vcfFile]['uri']

        # If the uri begins with 'file://', this needs to be removed
        finalUri = uri[6:] if uri.startswith('file://') else uri
        vcfName  = vcfFiles[vcfFile]['vcf_sample_name']
        if finalUri not in allVcfs: allVcfs.append(finalUri)
      mosaicSamples[sample]['vcf_file']        = finalUri
      mosaicSamples[sample]['vcf_sample_name'] = vcfName
      print('  Sample ', sample, ' appears as "', vcfName, '" in the header of vcf file: ', finalUri, sep = '')

  # If there are multiple VCF files, not all samples are in a single joint-called vcf. This case
  # is not yet handled:
  if len(allVcfs) > 1: fail('calypso_vcf_files.py: Multiple vcf files for the samples. Calypso requires a single multi sample vcf')

  # Return the updated mosaicSamples with the vcf information
  return mosaicSamples

# Check the reference genome of the vcf matches the project. This is crude check based on the length
# of chromosome 1
def checkReference(bcftools, reference, vcf):

  # Get the vcf header
  header = os.popen(bcf.getHeader(bcftools, vcf)).read()
  for line in header.split('\n'):
    if line.startswith('##contig=<ID=chr1,'):
      chr1Line = line.rstrip()
      break
    elif line.startswith('##contig=<ID=1,'):
      chr1Line = line.rstrip()
      break
    else: chr1Line = False

  # Check the length of chr1
  isCorrect    = False
  if 'length' not in chr1Line or not chr1Line: print('  Could not verify the reference genome for the vcf file')
  else:
    length = (line.rstrip().split('=')[3]).rstrip('>')
    if reference == 'GRCh38' and str(length).startswith('248956422'): isCorrect = True
    elif reference == 'GRCh37' and str(length).startswith('249250621'): isCorrect = True

  if isCorrect: print('  The vcf file matches the reference genome of the project (', str(reference), ')', sep = '')
  else:
    print('  The vcf file DOES NOT match the reference genome of the project (', str(reference), ')', sep = '')
    print('  Please verify the reference genome in Mosaic and ensure these match')
    exit(0)

# Determine the format of the chromosomes in the vcf file
def determineChromosomeFormat(bcftools, vcf):
  header = os.popen(bcf.getHeader(bcftools, vcf)).read()
  for line in header.split('\n'):
    if line.startswith('##contig'):
      chrId = (line.split('=')[2]).split(',')[0]
      if chrId.startswith('chr'): chrFormat = True
      else: chrFormat = False
      #else: fail('The first chromosome in the VCF file is "' + str(chrId) + '", not "1" or "chr1". This VCF file may not be sorted')
      break

  # Return the format
  return chrFormat

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
