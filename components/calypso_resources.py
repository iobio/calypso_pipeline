#!/usr/bin/python

from __future__ import print_function
from datetime import date
from os.path import exists
from sys import path

import os
import argparse
import json
import math
import glob
import importlib

# Check the supplied arguments
def checkResources(reference, dataDir, resourceFilename):
  resourceInfo = {}
  resourceInfo['path']   = dataDir + "GRCh" + str(reference) + "/"
  if not exists(resourceFilename): fail('The resource file "' + resourceFilename + '" does not exist')
  resourceInfo['json']   = resourceFilename

  # Return the info on resources
  return resourceInfo

# Read the json file describing the resources for the selected genome build, check the files exist, and store the versions.
def readResources(reference, resourceInfo):

  # Try and open the file
  try: resourceFile = open(resourceInfo['json'], 'r')
  except: fail('The file describing the resource files (' + str(resourceInfo['json']) + ') could not be found')

#  # Extract the json information
  try: resourceData = json.loads(resourceFile.read())
  except: fail('The json file (' + str(resourceFilename) + ') is not valid')

  # Store the resource data version
  try: resourceInfo['version'] = resourceData['version']
  except: fail('The resource json file (' + str(resourceInfo['json']) + ') does not include a version')

  # Store the date that this version was created
  try: resourceInfo['date'] = resourceData['date']
  except: fail('The resource json file (' + str(resourceInfo['json']) + ') does not include a date of creation')

  # Check that the resource json reference matches the selected reference
  try: resourceInfo['reference'] = resourceData['reference']
  except: fail('The resource json does not include a reference genome')
  if str(reference) != str(resourceData['reference']): fail("The selected reference (" + str(reference) + ") does not match the resource json reference (" + str(resourceInfo['reference']) + ")")

  # Get the resources
  try: resources = resourceData['resources']
  except: fail('The resources json does not include a list of resources')
  resourceInfo['resources'] = {}
  for resource in resources:
    resourceInfo['resources'][resource] = {}

    # Get required information for each resource
    try: resourceInfo['resources'][resource]['version'] = resources[resource]['version']
    except: fail('Version for resource "' + str(resource) + '" was not included in the resources json')

    # Get the associated resource file and attach the data directory to the file.
    try: resourceInfo['resources'][resource]['file'] = resources[resource]['file']
    except: fail("File for resource \"" + str(resource) + "\" was not included in the resources json")
    if resourceInfo['resources'][resource]['file']: resourceInfo["resources"][resource]["file"] = resourceInfo["path"] + resources[resource]["file"]

    # Annotation resources that are added using vcfanno need to be included in a toml file. The json
    # needs to include a Boolean to identify these
    try: resourceInfo['resources'][resource]['toml'] = resources[resource]['toml']
    except: fail('The resources json did not indicate if resource "' + str(resource) + '" should be included in the toml file')

    # If the resource is to be included in a toml file, instructions are required for vcfanno. The resouces COULD include
    # columns, fields, and names, but MUST include ops
    if resourceInfo['resources'][resource]['toml']:
      if 'fields' in resources[resource]: resourceInfo['resources'][resource]['fields'] = resources[resource]['fields']
      if 'columns' in resources[resource]: resourceInfo['resources'][resource]['columns'] = resources[resource]['columns']
      if 'names' in resources[resource]: resourceInfo['resources'][resource]['names'] = resources[resource]['names']
      try: resourceInfo['resources'][resource]['ops'] = resources[resource]['ops']
      except: fail('The resources json did not indicate the toml "ops" commands for resource "' + str(resource) + '"')

    # Check that the specified resource file exists
    if not exists(resourceInfo['resources'][resource]['file']): fail('Resource file ' + str(resourceInfo['resources'][resource]['file']) + ' does not exist')

    # Check if there are any postannotation commands for vcfanno
    if 'post_annotation' in resources[resource]:
      resourceInfo['resources'][resource]['post_annotation'] = {}
      for postAnnoField in resources[resource]['post_annotation']:
        if postAnnoField == 'name': resourceInfo['resources'][resource]['post_annotation']['name'] = resources[resource]['post_annotation']['name']
        elif postAnnoField == 'fields': resourceInfo['resources'][resource]['post_annotation']['fields'] = resources[resource]['post_annotation']['fields']
        elif postAnnoField == 'op': resourceInfo['resources'][resource]['post_annotation']['op'] = resources[resource]['post_annotation']['op']
        elif postAnnoField == 'type': resourceInfo['resources'][resource]['post_annotation']['type'] = resources[resource]['post_annotation']['type']
        else: fail('Unexpected post_annotation field in the resources json for resource "' + str(resource) + '"')

    # If the resource is vep, handle this separately
    if resource == 'vep':

      # VEP requires a cache and plugins directory to run. Get these directories and check they exist
      try: resourceInfo['resources']['vep']['cache'] = resources['vep']['cache']
      except: fail('The VEP resource description does not include the "cache"')
      try: resourceInfo['resources']['vep']['plugins'] = resources['vep']['plugins']
      except: fail('The VEP resource description does not include "plugins"')
    
      # Check that the supplied directories exist
      if not exists(resourceInfo['resources']['vep']['cache']): fail('VEP cache does not exist')
      if not exists(resourceInfo['resources']['vep']['plugins']): fail('VEP plugins directory does not exist')

  # Close the resources json
  resourceFile.close()

  # Return the updated resourceInfo
  return resourceInfo

# Output a summary file containing information on all of the resoources used to annotate the vcf
def calypsoSummary(workingDir, version, resourceInfo, reference):

  # Open an output summary file
  try: summaryFile = open(workingDir + 'calypso_' + str(date.today()) + '.txt', 'w')
  except: fail('Failed to open summary file')

  # Write relevant information to file
  print('### Output from Calypso pipeline ###', file = summaryFile)
  print(file = summaryFile)
  print('Calypso pipeline version: ', version, sep = '', file = summaryFile)
  print('Calypso resource version: ', resourceInfo['version'], sep = '', file = summaryFile)
  print('Reference:                ', reference, sep = '', file = summaryFile)
  print('Created on:               ', str(date.today()), sep = '', file = summaryFile)
  print(file = summaryFile)

  # Loop over all the used resources and output their versions
  for resource in resourceInfo['resources']:
    print(resource, file = summaryFile)
    print('  version: ', resourceInfo['resources'][resource]['version'], sep = '', file = summaryFile)
    print('  file:    ', resourceInfo['resources'][resource]['file'], sep = '', file = summaryFile)
    print(file = summaryFile)

  # Close the file
  summaryFile.close()

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
