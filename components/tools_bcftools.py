#!/usr/bin/python

from __future__ import print_function
from datetime import date
from os.path import exists
from sys import path

# Use bcftools query to return records in vcf format, but with only the selected INFO fields
def queryVcf(vcf, tags):
  command = 'bcftools query -f \'%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER'
  for tag in tags: command += '\\t%INFO/' + str(tag)
  command += '\\n\' ' + str(vcf)
  return command

# Use bcftools query to return the coords, alleles and selected info annotation fields
def query(vcf, tags):
  command = 'bcftools query -f \'%CHROM\\t%POS\\t%END\\t%REF\\t%ALT'
  for tag in tags: command += '\\t%INFO/' + str(tag)
  command += '\\n\' ' + str(vcf)
  return command

# Use bcftools query to return the coords, alleles and "worst" consequence CSQ fields
def splitvepQuery(vcf, tags):
  command = 'bcftools +split-vep -s worst -c CSQ -f \'%CHROM\\t%POS\\t%END\\t%REF\\t%ALT'
  for tag in tags: command += '\\t%INFO/' + str(tag)
  command += '\\n\' ' + str(vcf)
  return command

# Use bcftools query to return the coords, alleles and selected format (genotype) annotation fields
def queryFormat(vcf, tag):
  command = 'bcftools query -f \'%CHROM\\t%POS\\t%END\\t%REF\\t%ALT[\\t%' + str(tag) + ']'
  command += '\\n\' ' + str(vcf)
  return command

# Generate intersections and complements between two vcf files
def isec(oldVcf, newVcf, cvDir):
  command = 'bcftools isec -O z -p ' + str(cvDir) + ' ' + str(oldVcf) + ' ' + str(newVcf)
  return command

# Generate the intersection between two vcf files
def isecInt(vcfA, vcfB, outVcf):
  command = 'bcftools isec -n =2 -w 1 -O z -o ' + str(outVcf) + ' ' + str(vcfA) + ' ' + str(vcfB)
  return command

# Generate the intersection between two vcf files and stream the output
def isecIntStream(vcfA, vcfB):
  command = 'bcftools isec -n =2 -w 1 ' + str(vcfA) + ' ' + str(vcfB)
  return command

# Generate the intersection between two vcf files, but stream the output with no header
def isecIntNoHeader(vcfA, vcfB):
  command = 'bcftools isec -n =2 -w 1 -H ' + str(vcfA) + ' ' + str(vcfB)
  return command

# Use bcftools view to get a vcf header
def getHeader(inVcf, outVcf):
  command = 'bcftools view -h -O v -o ' + str(outVcf) + ' ' + str(inVcf)
  return command

# Compress a vcf file
def compress(vcf, outType):
  return 'bcftools view -O ' + str(outType) + ' -o ' + str(vcf) + '.gz ' + str(vcf)

# Index a file
def index(vcf):
  return 'bcftools index -t ' + str(vcf)

# Merge a list of vcf files
def merge(inVcfs, outType, outVcf):
  command = 'bcftools merge -O ' + str(outType) + ' -o ' + str(outVcf)
  for inFile in inFiles: command += ' ' + str(inFile)
  return command

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
