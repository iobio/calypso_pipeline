#!/usr/bin/python

from datetime import date
from os.path import exists
from sys import path

import os

# Determine if the defined executable can be used
def isExecutable(bcftoolsExe):
  return os.access(bcftoolsExe, os.X_OK)

# Get the bcftools version
def version(bcftools):
  command = bcftools + ' -v'
  return command

#####
##### bcftools view commands
#####

# Use bcftools view to generate a vcf header as a vcf file
def createHeaderVcf(bcftools, inVcf, outVcf):
  command = bcftools + ' view -h -O v -o ' + str(outVcf) + ' ' + str(inVcf)
  return command

# Use bcftools view to return a vcf header
def getHeader(bcftools, vcf):
  command = bcftools + ' view -h ' + str(vcf)
  return command

# Compress a vcf file
def compress(bcftools, inVcf, outVcf, outType):
  return bcftools + ' view -O ' + str(outType) + ' -o ' + str(outVcf) + ' ' + str(inVcf)

# Compress a vcf file and create the index
def compressAndIndex(bcftools, inVcf, outVcf, outType):
  return bcftools + ' view -O ' + str(outType) + ' -o ' + str(outVcf) + ' --write-index ' + str(inVcf)

#####
##### bcftools query commands
#####

# Use bcftools query to return records in vcf format, but with only the selected INFO fields
def queryVcf(bcftools, vcf, tags):
  command = bcftools + ' query -f \'%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER'
  for tag in tags: command += '\\t%INFO/' + str(tag)
  command += '\\n\' ' + str(vcf)
  return command

# Use bcftools query to return the coords, alleles and selected info annotation fields
def query(bcftools, vcf, tags):
  command = bcftools + ' query -f \'%CHROM\\t%POS\\t%END\\t%REF\\t%ALT'
  for tag in tags: command += '\\t%INFO/' + str(tag)
  command += '\\n\' ' + str(vcf)
  return command

# Use bcftools query to return the coords, alleles and selected format (genotype) annotation fields
def queryFormat(bcftools, vcf, tag):
  command = bcftools + ' query -f \'%CHROM\\t%POS\\t%END\\t%REF\\t%ALT[\\t%' + str(tag) + ']'
  command += '\\n\' ' + str(vcf)
  return command

#####
##### bcftools isec commands
#####

# Generate intersections and complements between two vcf files
def isec(bcftools, oldVcf, newVcf, cvDir):
  command = bcftools + ' isec -O z -p ' + str(cvDir) + ' ' + str(oldVcf) + ' ' + str(newVcf)
  return command

# Generate the intersection between two vcf files
def isecInt(bcftools, vcfA, vcfB, outVcf):
  command = bcftools + ' isec -n =2 -w 1 -O z -o ' + str(outVcf) + ' ' + str(vcfA) + ' ' + str(vcfB)
  return command

# Generate the intersection between two vcf files and stream the output
def isecIntStream(bcftools, vcfA, vcfB):
  command = bcftools + ' isec -n =2 -w 1 ' + str(vcfA) + ' ' + str(vcfB)
  return command

# Generate the intersection between two vcf files, but stream the output with no header
def isecIntNoHeader(bcftools, vcfA, vcfB):
  command = bcftools + ' isec -n =2 -w 1 -H ' + str(vcfA) + ' ' + str(vcfB)
  return command

#####
##### bcftools index commands
#####

# Index a file
def index(bcftools, vcf):
  return bcftools + ' index -t ' + str(vcf)

#####
##### bcftools merge commands
#####

# Merge a list of vcf files
def merge(bcftools, inVcfs, outType, outVcf):
  command = bcftools + ' merge -O ' + str(outType) + ' -o ' + str(outVcf)
  command += ' --write-index'
  for inVcf in inVcfs: command += ' ' + str(inVcf)
  return command

#####
##### bcftools annotate commands
#####

# Extract annotations from a vcf file
def extractAnnotations(bcftools, tags, addedAnnotations, header, renameFile, inVcf, outVcf):
  command = bcftools + ' annotate -x '

  # Add the annotations to keep
  for i, tag in enumerate(tags):
    if i == 0: command += '^INFO/' + str(tag)
    else: command += ',^INFO/' + str(tag)

  # If a header line is provided, add it
  if header: command += ' --header-lines ' + str(header)

  # If a rename file is included, add this to the command
  if renameFile: command += ' --rename-annots ' + str(renameFile)

  # If a list of annotations to add is included, add them
  if addedAnnotations:
    command += ' -m '
    for annotation in addedAnnotations: command += str(annotation) + ','
    command = command.rstrip(',')

  # If an output file is included, output this as a vcf.gz file, otherwise assume this command will be
  # piped to another command
  if outVcf: command += ' -O z -o ' + str(outVcf) + ' --write-index'

  # Finally add the input vcf file
  command += ' ' + str(inVcf)

  return command

# Annotate a vcf with annotations from a different vcf
def annotateWithVcf(bcftools, tags, annotateVcf, inVcf, outVcf):
  command = bcftools + ' annotate -a ' + str(annotateVcf) + ' -c '

  # Add the annotations to keep
  for i, tag in enumerate(tags):
    if i == 0: command += 'INFO/' + str(tag)
    else: command += ',INFO/' + str(tag)

  # If an output file is included, output this as a vcf.gz file, otherwise assume this command will be
  # piped to another command
  if outVcf: command += ' -O z -o ' + str(outVcf) + ' --write-index'

  # Finally add the input vcf file
  command += ' ' + str(inVcf)

  return command

#####
##### bcftools split-vep commands
#####

# Use bcftools query to return the coords, alleles and "worst" consequence CSQ fields
def splitvepQuery(bcftools, vcf, tags):
  command = bcftools + ' +split-vep -s worst -c CSQ -f \'%CHROM\\t%POS\\t%END\\t%REF\\t%ALT'
  for tag in tags: command += '\\t%INFO/' + str(tag)
  command += '\\n\' ' + str(vcf)
  return command


# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
