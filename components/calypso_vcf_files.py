from __future__ import print_function

import os
import json
import tools_bcftools as bcf

# Get the vcf files for all the mosaic samples
def getVcfFiles(mosaicConfig, api_sf, projectId, mosaicSamples):
  allVcfs = []

  # Loop over all the samples
  for sample in mosaicSamples:
    sampleId = mosaicSamples[sample]['id']
    vcfFiles = api_sf.getSampleVcfs(mosaicConfig, projectId, sampleId)

    # If there are no vcf files for this sample output a warning, but continue, unless this is the proband
    if len(vcfFiles) == 0:
      if mosaicSamples[sample]['relation'] == 'Proband': fail('calypso_vcf_files: The proband (' + sample + ') is not associated with any vcf files. The proband must appear in a vcf file')
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

# Determine the format of the chromosomes in the vcf file
def determineChromosomeFormat(bcftools, vcf):
  header = os.popen(bcf.getHeader(bcftools, vcf)).read()
  for line in header.split('\n'):
    if line.startswith('##contig'):
      chrId = (line.split('=')[2]).split(',')[0]
      if chrId == 'chr1': chrFormat = True
      elif chrId == '1': chrFormat = False
      else: fail('The first chromosome in the VCF file is "' + str(chrId) + '", not "1" or "chr1". This VCF file may not be sorted')
      break

  # Return the format
  return chrFormat

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
