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
import tools_bcftools as bcf

import subprocess
from subprocess import Popen

# Check the supplied arguments
def checkResources(reference, dataDir, toolsDir, resourceFilename):
  resourceInfo = {}
  resourceInfo['path'] = dataDir + str(reference) + '/'

  # Store the directory where tools reside if defined
  if toolsDir: resourceInfo['toolsPath'] = str(toolsDir) if toolsDir.endswith('/') else str(toolsDir) + str('/')
  else: resourceInfo['toolsPath'] = False

  if not exists(resourceFilename): fail('The resource file "' + resourceFilename + '" does not exist')
  resourceInfo['json']   = resourceFilename

  # Return the info on resources
  return resourceInfo

# Read the json file describing the resources for the selected genome build, check the files exist, and store the versions.
def readResources(reference, rootPath, resourceInfo):

  # Try and open the file
  try: resourceFile = open(resourceInfo['json'], 'r')
  except: fail('The file describing the resource files (' + str(resourceInfo['json']) + ') could not be found')

  # Extract the json information
  try: resourceData = json.loads(resourceFile.read())
  except: fail('The json file (' + str(resourceInfo['json']) + ') is not valid')

  # Check that this is a resources file and not a Mosaic definition file
  try: resourceInfo['resourceType'] = resourceData['type']
  except: fail('The resource json file (' + str(resourceInfo['json']) + ') does not include a type (calypso_resource, or mosaic_resource)')
  if str(resourceInfo['resourceType']) != 'calypso_resource': fail('The specified resource file (' + str(resourceInfo['json']) + ') is not specified as a "calypso_resource" in the "type" field')

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

  # Check if a path to a tools directory is provided. If the path was set on the command line, use the value provided on the command line
  if not resourceInfo['toolsPath'] and 'tools_path' in resourceData: 
    resourceInfo['toolsPath'] = rootPath + '/' + resourceData['tools_path']
    if not resourceInfo['toolsPath'].endswith('/'): resourceInfo['toolsPath'] += '/'

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

      # If the vep cache / plugins directories include "TOOLS", this should be replaced with the path to the
      # tools
      if 'TOOLS' in resourceInfo['resources']['vep']['cache']:
        if not resourceInfo['toolsPath']: continue #fail('The resources json includes "TOOLS" in the vep cache path, which requires a tools path be defined on the command line')
        resourceInfo['resources']['vep']['cache']   = resourceInfo['resources']['vep']['cache'].replace('TOOLS/', resourceInfo['toolsPath'])
        resourceInfo['resources']['vep']['plugins'] = resourceInfo['resources']['vep']['plugins'].replace('TOOLS/', resourceInfo['toolsPath'])
        if resourceInfo['resources']['vep']['plugins'].endswith('/'): resourceInfo['resources']['vep']['plugins'].rstrip('/')

      # Check that the supplied directories exist
      if not exists(resourceInfo['resources']['vep']['cache']): fail('VEP cache does not exist (' + str(resourceInfo['resources']['vep']['cache']) + ')')
      if not exists(resourceInfo['resources']['vep']['plugins']): fail('VEP plugins directory does not exist (' + str(resources['resources']['vep']['plugins']) + ')')

      # Check if any plugins are to be used. These would be stored in "active_plugins"
      if 'active_plugins' in resources['vep']:
        resourceInfo['resources']['vep']['active_plugins'] = {}
        for plugin in resources['vep']['active_plugins']:
          resourceInfo['resources']['vep']['active_plugins'][plugin] = {}

          # Get the files for the plugin
          pluginFiles = ''
          if 'files' in resources['vep']['active_plugins'][plugin]:
            for i, pluginFile in enumerate(resources['vep']['active_plugins'][plugin]['files']):
              if 'PLUGIN' in pluginFile: pluginFile = str(pluginFile.replace('PLUGIN', str(resourceInfo['resources']['vep']['plugins'])))
              if 'CACHE' in pluginFile: pluginFile = str(pluginFile.replace('CACHE', str(resourceInfo['resources']['vep']['cache'])))

              # The resource for the VEP plugin may be stored in the Calypso data store. If so, get the name of the resource and determine
              # where the data lives
              if 'RESOURCE' in pluginFile:

                # Get the name of the resource
                if ':' not in pluginFile: fail('VEP plugin "' + plugin + '" has a required file associated with a resource. This must be specified as "RESOURCE:resource_name"')
                resourceName = pluginFile.split(':')[1]
                if resourceName not in resourceInfo['resources']: fail('VEP plugin "' + plugin + '" has a required file associated with resource "' + resourceName + '" which does not exist')
                if 'file' not in resourceInfo['resources'][resourceName]: fail('VEP plugin "' + plugin + '" has a required file associated with resource "' + resourceName + '" which does not have any associated files')
                pluginFile = pluginFile.split('RESOURCE')[0] + resourceInfo['resources'][resourceName]['file']
              if i == 0: pluginFiles += pluginFile
              else: pluginFiles += ',' + pluginFile
          if pluginFiles != '': resourceInfo['resources']['vep']['active_plugins'][plugin]['files'] = pluginFiles

          # If a path needs to be updated, it will included in an 'export' command
          if 'export' in resources['vep']['active_plugins'][plugin]:
            if 'export' not in resourceInfo['resources']['vep']: resourceInfo['resources']['vep']['export'] = []
            for newPath in resources['vep']['active_plugins'][plugin]['export']:
              if 'CACHE' in newPath: newPath = newPath.replace('CACHE', str(resourceInfo['resources']['vep']['cache']))
              if 'PLUGIN' in newPath: newPath = newPath.replace('PLUGIN', str(resourceInfo['resources']['vep']['plugins']))
            resourceInfo['resources']['vep']['export'].append('export ' + newPath)

      # Check for additional options to be included in the VEP command
      if 'options' in resources['vep']:
        if type(resources['vep']['options']) != list: fail('VEP options need to be included in the resources file as a list')
        resourceInfo['resources']['vep']['options'] = resources['vep']['options']

      # The fields that are to be included in the VEP output needs to be defined
      if 'fields' not in resources['vep']: fail('VEP fields need to be defined')
      if type(resources['vep']['fields']) != list: fail('The CSQ annotations in the "fields" section of the VEP resource description must be a list')
      resourceInfo['resources']['vep']['fields'] = resources['vep']['fields']

      # The CSQ fields that need to be extracted for use by Slivar need to be defined in the resource file
      if 'slivar' not in resources['vep']: fail('VEP CSQ annotations that need to be extracted for use by slivar need to be defined')
      if type(resources['vep']['slivar']) != list: fail('The CSQ annotations in the "slivar" section of the VEP resource description must be a list')
      isFail = False
      missingAnnotations = []
      for annotation in resources['vep']['slivar']:

        # The annotation can include a definition of type. Remove this and check the annotation name is included in the "fields" VEP
        # argument. Any annotation that is to be extracted for slivar has to exist in the VEP output
        name = annotation.split(':')[0] if ':' in annotation else annotation
        if name not in resourceInfo['resources']['vep']['fields']:
          missingAnnotations.append(name)
          isFail = True
      if isFail: fail('The following VEP annotations are listed as required by Slivar, but are not included in the fields to be output by VEP:\n  ' + ', '.join(missingAnnotations))
      resourceInfo['resources']['vep']['slivar'] = resources['vep']['slivar']

  # Close the resources json
  resourceFile.close()

  # Return the updated resourceInfo
  return resourceInfo

# Define all of the tools to be used in Calypso, check they exist and output their versions
def calypsoTools(resourceInfo):

  # Define the tools. If no path to the tools is provided, assume the tools are available in the system and don't need a path.
  # The preference is to use the parh to the Calypso tools directory
  resourceInfo['tools'] = {}
  if resourceInfo['toolsPath']:
    resourceInfo['tools']['bcftools'] = str(resourceInfo['toolsPath']) + 'bcftools/bcftools'
    resourceInfo['tools']['vcfanno']  = str(resourceInfo['toolsPath']) + 'vcfanno'
    resourceInfo['tools']['slivar']   = str(resourceInfo['toolsPath']) + 'slivar'
    resourceInfo['tools']['vep']      = str(resourceInfo['toolsPath']) + 'ensembl-vep/vep'
  else:
    resourceInfo['tools']['bcftools'] = 'bcftools'
    resourceInfo['tools']['vcfanno']  = 'vcfanno'
    resourceInfo['tools']['slivar']   = 'slivar'
    resourceInfo['tools']['vep']      = 'vep'
  print('  Using the following tools:')

  # Bcftools
  bcftoolsVersion = os.popen(bcf.version(resourceInfo['tools']['bcftools'])).readlines()[0].rstrip()
  print('    ', bcftoolsVersion, sep = '')

  # vcfanno
  proc = Popen(resourceInfo['tools']['vcfanno'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  print('    ', proc.stderr.readlines()[2].decode().rstrip(), sep = '')

  # slivar
  proc = Popen(resourceInfo['tools']['slivar'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  print('    ', proc.stderr.readlines()[0].decode().rstrip()[2:], sep = '')

  # vep
  vepInfo = os.popen(resourceInfo['tools']['vep']).readlines()
  print('    vep versions:')
  for i in range(5, 9, 1): print('    ', vepInfo[i].rstrip(), sep = '')

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
