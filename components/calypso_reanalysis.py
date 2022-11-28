#!/usr/bin/python

from __future__ import print_function
from datetime import date
from os.path import exists
from sys import path

import importlib
import json
import os
import tools_bcftools as bcf

# Compare the current and previous resources and identify any differences
def compare(previousResourceInfo, resourceInfo):
  addedResources   = []
  removedResources = list(previousResourceInfo['resources'].keys())
  updatedResources = {}

  # Loop over the resources in the current resource json and compare to the previous json
  for resource in resourceInfo['resources']:
    if resource not in previousResourceInfo['resources']: addedResources.append(resource)
    else:
      if str(resourceInfo['resources'][resource]['version']) != str(previousResourceInfo['resources'][resource]['version']):
        updatedResources[resource] = {'previousVersion': previousResourceInfo['resources'][resource]['version'], 'currentVersion': resourceInfo['resources'][resource]['version']}
      removedResources.remove(resource)

  # Write out a summary of changes
  if len(removedResources) > 0:
    print()
    print('The following resources have been removed from Calypso so will not be updated:')
    for resource in removedResources: print('  ', resource)
  if len(addedResources) > 0:
    print()
    print('The following resources have been added to Calypso, so this will be the first inclusion of these annotations:')
    for resource in addedResources: print('  ', resource)
  if len(updatedResources) > 0:
    print()
    print('The following resources have been updated:')
    for resource in updatedResources: print('  ', resource, ' was updated from ', updatedResources[resource]['previousVersion'], ' to ', updatedResources[resource]['currentVersion'], sep = '')

  # Return the updated resources
  return updatedResources

# Determine which project variants have undergone notable changes
def clinvarUpdates(cvDates, compDir, bashFile):

  # Define the name of the output file
  outFile = 'clinvar_updates_' + str(cvDates) + '.tsv'

  # Define the clinvar difference files to use for intersecting with the project vcf file
  clinvar = {}
  clinvar['benign_conflicting.vcf.gz'] = 'benign_conflicting'
  clinvar['benign_path.vcf.gz'] = 'benign_pathogenic'
  clinvar['benign_vus.vcf.gz'] = 'banign_vus'
  clinvar['conflicting_path.vcf.gz'] = 'conflicting_pathogenic'
  clinvar['conflicting_vus.vcf.gz'] = 'conflicting_vus'
  clinvar['path_conflicting.vcf.gz'] = 'pathogenic_conflicting'
  clinvar['path_vus.vcf.gz'] = 'pathogenic_vus'
  clinvar['vus_conflicting.vcf.gz'] = 'vus_conflicting'
  clinvar['vus_path.vcf.gz'] = 'vus_pathogenic'
  clinvar['unique_current_conflicting.vcf.gz'] = 'new_conflicting'
  clinvar['unique_current_path.vcf.gz'] = 'new_pathogenic'
  clinvar['unique_current_vus.vcf.gz'] = 'new_vus'
  clinvar['unique_previous_conflicting.vcf.gz'] = 'removed_conflicting'
  clinvar['unique_previous_path.vcf.gz'] = 'removed_pathogenic'
  clinvar['unique_previous_vus.vcf.gz'] = 'removed_vus'
  print(file = bashFile)
  print('# ClinVar updates', file = bashFile)
  print('echo -n "Searching for ClinVar updates..."', file = bashFile)
  for cvFile in clinvar:
    print(bcf.isecIntStream(str(compDir) + '/' + cvFile, '$FINALVCF'), ' \\', file = bashFile)
    print(' | bcftools query -f \'%CHROM\\t%POS\\t%END\\t%REF\\t%ALT\\t', clinvar[cvFile], '\\n \\', sep = '', file = bashFile)
    print(' >> ', outFile, sep = '', file = bashFile)
  print('echo "complete"', file = bashFile)

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
