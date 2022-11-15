#!/usr/bin/python

from __future__ import print_function
from datetime import date
from os.path import exists
from sys import path

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

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
