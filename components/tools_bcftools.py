#!/usr/bin/python

from __future__ import print_function
from datetime import date
from os.path import exists
from sys import path

# Use bcftools query to return records in vcf format, but with only the selected INFO fields
def queryVcf(vcf, tags):

  # Build the bcftools query command beginning with the base command
  command = 'bcftools query -f \'%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER'

  # Add all of the custom INFO tags
  for tag in tags: command += '\\t%INFO/' + str(tag)

  # Complete the command with a newline and the vcf file
  command += '\\n\' ' + str(vcf)

  # Return the complete command
  return command

# Use bcftools query to return the coords, alleles and selected info annotation fields
def query(vcf, tags):

  # Build the bcftools query command beginning with the base command
  command = 'bcftools query -f \'%CHROM\\t%POS\\t%END\\t%REF\\t%ALT'

  # Add all of the custom INFO tags
  for tag in tags: command += '\\t%INFO/' + str(tag)

  # Complete the command with a newline and the vcf file
  command += '\\n\' ' + str(vcf)

  # Return the complete command
  return command

# Use bcftools query to return the coords, alleles and "worst" consequence CSQ fields
def splitvepQuery(vcf, tags):

  # Build the bcftools query command beginning with the base command
  command = 'bcftools +split-vep -s worst -c CSQ -f \'%CHROM\\t%POS\\t%END\\t%REF\\t%ALT'

  # Add all of the custom INFO tags
  for tag in tags: command += '\\t%INFO/' + str(tag)

  # Complete the command with a newline and the vcf file
  command += '\\n\' ' + str(vcf)

  # Return the complete command
  return command

# Use bcftools query to return the coords, alleles and selected format (genotype) annotation fields
def queryFormat(vcf, tag):

  # Build the bcftools query command beginning with the base command
  command = 'bcftools query -f \'%CHROM\\t%POS\\t%END\\t%REF\\t%ALT[\\t%' + str(tag) + ']'

  # Complete the command with a newline and the vcf file
  command += '\\n\' ' + str(vcf)

  # Return the complete command
  return command

# Generate intersections and complements between two vcf files
def isec(oldVcf, newVcf, cvDir):

  # Build the bcftools command
  command = 'bcftools isec -O z -p ' + str(cvDir) + ' ' + str(oldVcf) + ' ' + str(newVcf)

  # Return the command
  return command

# Generate the intersection between two vcf files
def isecInt(vcfA, vcfB, outVcf):

  # Build the bcftools command
  command = 'bcftools isec -n =2 -w 1 -O z -o ' + str(outVcf) + ' ' + str(vcfA) + ' ' + str(vcfB)

  # Return the command
  return command

# Use bcftools view to get a vcf header
def getHeader(vcf):

  # Build the bcftools view command
  command = 'bcftools view -h ' + str(vcf)

  # Return the complete command
  return command

# Compress a vcf file
def compress(vcf, outputType):
  return 'bcftools view -O ' + str(outputType) + ' -o ' + str(vcf) + '.gz ' + str(vcf)

# Index a file
def index(vcf):
  return 'bcftools index -t ' + str(vcf)

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
