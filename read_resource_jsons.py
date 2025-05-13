import json
import os
import subprocess
import sys

from os.path import exists
from sys import path

# Read the json file describing the resources for the selected genome build, check the files exist, and store the versions.
def read_resources(reference, root_path, resource_info, use_vep):

  # Try and open the file
  try:
    resource_file = open(resource_info['json'], 'r')
  except:
    fail('The file describing the resource files (' + str(resource_info['json']) + ') could not be found')

  # Extract the json information
  try:
    resource_data = json.loads(resource_file.read())
  except:
    fail('The json file (' + str(resource_info['json']) + ') is not valid')

  # Check that this is a resources file and not a Mosaic definition file
  try:
    resource_info['resourceType'] = resource_data['type']
  except:
    fail('The resource json file (' + str(resource_info['json']) + ') does not include a type (calypso_resource, or mosaic_resource)')
  if str(resource_info['resourceType']) != 'calypso_resource':
    fail('The specified resource file (' + str(resource_info['json']) + ') is not specified as a "calypso_resource" in the "type" field')

  # Store the resource data version
  try:
    resource_info['version'] = resource_data['version']
  except:
    fail('The resource json file (' + str(resource_info['json']) + ') does not include a version')

  # Store the date that this version was created
  try:
    resource_info['date'] = resource_data['date']
  except:
    fail('The resource json file (' + str(resource_info['json']) + ') does not include a date of creation')

  # Check that the resource json reference matches the selected reference
  try:
    resource_info['reference'] = resource_data['reference']
  except:
    fail('The resource json does not include a reference genome')
  if str(reference) != str(resource_data['reference']):
    fail("The selected reference (" + str(reference) + ") does not match the resource json reference (" + str(resource_info['reference']) + ")")

  # Check if an SV json is specified
  try:
    resource_info['sv_json'] = resource_data['SV json']
  except:
    resource_info['sv_json'] = False

  # Check if a path to a tools directory is provided. If the path was set on the command line, use the value provided on the command line
  if not resource_info['toolsPath'] and 'tools_path' in resource_data:
    resource_info['toolsPath'] = root_path + '/' + resource_data['tools_path']
    if not resource_info['toolsPath'].endswith('/'):
      resource_info['toolsPath'] += '/'

  # Get the resources
  try:
    resources = resource_data['resources']
  except:
    fail('The resources json does not include a list of resources')
  resource_info['resources'] = {}
  for resource in resources:
    resource_info['resources'][resource] = {}

    # Get required information for each resource
    try:
      resource_info['resources'][resource]['version'] = resources[resource]['version']
    except:
      fail('Version for resource "' + str(resource) + '" was not included in the resources json')

    # Get the associated resource file and attach the data directory to the file.
    try:
      resource_info['resources'][resource]['file'] = resources[resource]['file']
    except:
      fail('File for resource "' + str(resource) + '" was not included in the resources json')

    # Annotation resources that are added using vcfanno need to be included in a toml file. The json
    # needs to include a Boolean to identify these
    try:
      resource_info['resources'][resource]['toml'] = resources[resource]['toml']
    except:
      fail('The resources json did not indicate if resource "' + str(resource) + '" should be included in the toml file')

    # If the resource is to be included in a toml file, instructions are required for vcfanno. The resouces COULD include
    # columns, fields, and names, but MUST include ops
    if resource_info['resources'][resource]['toml']:
      if 'fields' in resources[resource]:
        resource_info['resources'][resource]['fields'] = resources[resource]['fields']
      if 'columns' in resources[resource]:
        resource_info['resources'][resource]['columns'] = resources[resource]['columns']
      if 'names' in resources[resource]:
        resource_info['resources'][resource]['names'] = resources[resource]['names']
      try:
        resource_info['resources'][resource]['ops'] = resources[resource]['ops']
      except:
        fail('The resources json did not indicate the toml "ops" commands for resource "' + str(resource) + '"')

    # Check that the specified resource file exists
    if resource_info['resources'][resource]['file']:
      full_file_path = resource_info['path'] + resource_info['resources'][resource]['file']
      if not exists(full_file_path):
        fail('Resource file ' + str(resource_info['resources'][resource]['file']) + ' does not exist')

    # Check if there are any postannotation commands for vcfanno
    if 'post_annotation' in resources[resource]:
      resource_info['resources'][resource]['post_annotation'] = {}
      for postAnnoField in resources[resource]['post_annotation']:
        if postAnnoField == 'name':
          resource_info['resources'][resource]['post_annotation']['name'] = resources[resource]['post_annotation']['name']
        elif postAnnoField == 'fields':
          resource_info['resources'][resource]['post_annotation']['fields'] = resources[resource]['post_annotation']['fields']
        elif postAnnoField == 'op':
          resource_info['resources'][resource]['post_annotation']['op'] = resources[resource]['post_annotation']['op']
        elif postAnnoField == 'type':
          resource_info['resources'][resource]['post_annotation']['type'] = resources[resource]['post_annotation']['type']
        else:
          fail('Unexpected post_annotation field in the resources json for resource "' + str(resource) + '"')

    # If the resource is vep, handle this separately
    if resource == 'vep':
      if use_vep:
        resource_info['resources']['vep']['ignore'] = False
      else:
        resource_info['resources']['vep']['ignore'] = True

      # VEP requires a cache and plugins directory to run. Get these directories and check they exist
      try:
        resource_info['resources']['vep']['cache'] = resources['vep']['cache']
      except:
        fail('The VEP resource description does not include the "cache"')
      try:
        resource_info['resources']['vep']['plugins'] = resources['vep']['plugins']
      except:
        fail('The VEP resource description does not include "plugins"')

      # If the vep cache / plugins directories include "TOOLS", this should be replaced with the path to the
      # tools
      if 'TOOLS' in resource_info['resources']['vep']['cache']:
        if not resource_info['toolsPath']:
          fail('The resources json includes "TOOLS" in the vep cache path, which requires a tools path be defined on the command line')
        resource_info['resources']['vep']['cache'] = resource_info['resources']['vep']['cache'].replace('TOOLS/', resource_info['toolsPath'])
        resource_info['resources']['vep']['plugins'] = resource_info['resources']['vep']['plugins'].replace('TOOLS/', resource_info['toolsPath'])
        if resource_info['resources']['vep']['plugins'].endswith('/'):
          resource_info['resources']['vep']['plugins'].rstrip('/')

      # Check that the supplied directories exist
      if not exists(resource_info['resources']['vep']['cache']):
        fail('VEP cache does not exist (' + str(resource_info['resources']['vep']['cache']) + ')')
      if not exists(resource_info['resources']['vep']['plugins']):
        fail('VEP plugins directory does not exist (' + str(resources['resources']['vep']['plugins']) + ')')

      # Check if any plugins are to be used. These would be stored in "active_plugins"
      if 'active_plugins' in resources['vep']:
        resource_info['resources']['vep']['active_plugins'] = {}
        for plugin in resources['vep']['active_plugins']:
          resource_info['resources']['vep']['active_plugins'][plugin] = {}

          # Get the files for the plugin
          plugin_files = ''
          if 'files' in resources['vep']['active_plugins'][plugin]:
            for i, plugin_file in enumerate(resources['vep']['active_plugins'][plugin]['files']):
              if 'PLUGIN' in plugin_file:
                plugin_file = str(plugin_file.replace('PLUGIN', str(resource_info['resources']['vep']['plugins'])))
              if 'CACHE' in plugin_file:
                plugin_file = str(plugin_file.replace('CACHE', str(resource_info['resources']['vep']['cache'])))

              # The resource for the VEP plugin may be stored in the Calypso data store. If so, get the name of the resource and determine
              # where the data lives
              if 'RESOURCE' in plugin_file:

                # Get the name of the resource
                if ':' not in plugin_file:
                  fail('VEP plugin "' + plugin + '" has a required file associated with a resource. This must be specified as "RESOURCE:resource_name"')
                resourceName = plugin_file.split(':')[1]
                if resource_name not in resource_info['resources']:
                  fail('VEP plugin "' + plugin + '" has a required file associated with resource "' + resource_name + '" which does not exist')
                if 'file' not in resource_info['resources'][resource_name]:
                  fail('VEP plugin "' + plugin + '" has a required file associated with resource "' + resource_name + '" which does not have any associated files')
                plugin_file = plugin_file.split('RESOURCE')[0] + resource_info['resources'][resource_name]['file']
              if i == 0:
                plugin_files += plugin_file
              else:
                plugin_files += ',' + plugin_file
          if plugin_files != '':
            resource_info['resources']['vep']['active_plugins'][plugin]['files'] = plugin_files

          # If a path needs to be updated, it will included in an 'export' command
          if 'export' in resources['vep']['active_plugins'][plugin]:
            if 'export' not in resource_info['resources']['vep']:
              resource_info['resources']['vep']['export'] = []
            for new_path in resources['vep']['active_plugins'][plugin]['export']:
              if 'CACHE' in new_path:
                new_path = new_path.replace('CACHE', str(resource_info['resources']['vep']['cache']))
              if 'PLUGIN' in new_path:
                new_path = new_path.replace('PLUGIN', str(resource_info['resources']['vep']['plugins']))
            resource_info['resources']['vep']['export'].append('export ' + new_path)

      # Check for additional options to be included in the VEP command
      if 'options' in resources['vep']:
        if type(resources['vep']['options']) != list:
          fail('VEP options need to be included in the resources file as a list')
        resource_info['resources']['vep']['options'] = resources['vep']['options']

      # The fields that are to be included in the VEP output needs to be defined
      if 'fields' not in resources['vep']:
        fail('VEP fields need to be defined')
      if type(resources['vep']['fields']) != list:
        fail('The CSQ annotations in the "fields" section of the VEP resource description must be a list')
      resource_info['resources']['vep']['fields'] = resources['vep']['fields']

      # The CSQ fields that need to be extracted for use by Slivar need to be defined in the resource file
      if 'slivar' not in resources['vep']:
        fail('VEP CSQ annotations that need to be extracted for use by slivar need to be defined')
      if type(resources['vep']['slivar']) != list:
        fail('The CSQ annotations in the "slivar" section of the VEP resource description must be a list')
      is_fail = False
      missing_annotations = []
      for annotation in resources['vep']['slivar']:

        # The annotation can include a definition of type. Remove this and check the annotation name is included in the "fields" VEP
        # argument. Any annotation that is to be extracted for slivar has to exist in the VEP output
        name = annotation.split(':')[0] if ':' in annotation else annotation
        if name not in resource_info['resources']['vep']['fields']:
          missing_annotations.append(name)
          is_fail = True
      if is_fail:
        fail('The following VEP annotations are listed as required by Slivar, but are not included in the fields to be output by VEP:\n  ' + ', '.join(missing_annotations))
      resource_info['resources']['vep']['slivar'] = resources['vep']['slivar']

  # Close the resources json
  resource_file.close()

  # Return the updated resource_info
  return resource_info

######
###### Following are routines for parsing the mosaic json
######

def read_mosaic_json(filename, reference):
  mosaic_info = {}
  mosaic_info['resources'] = {}

  # Try and open the file
  try:
    mosaic_file = open(filename, 'r')
  except:
    fail('The file describing Mosaic related information (' + str(filename) + ') could not be found')

  # Extract the json information
  try:
    mosaic_data = json.loads(mosaic_file.read())
  except:
    fail('The json file (' + str(filename) + ') is not valid')

  # Check that this is a Mosaic resources file and not a Calypso resources file
  try:
    mosaic_info['resourceType'] = mosaic_data['type']
  except:
    fail('The resource json file (' + str(mosaic_info['json']) + ') does not include a type (calypso_resource, or mosaic_resource)')
  if str(mosaic_info['resourceType']) != 'mosaic_resource':
    fail('The specified resource file (' + str(mosaic_info['json']) + ') is not specified as a "mosaic_resource" in the "type" field')

  # Store the data version
  try:
    mosaic_info['version'] = mosaic_data['version']
  except:
    fail('The Mosaic json (' + str(filename) + ') does not include a version')

  # Store the date the version was created
  try:
    mosaic_info['date'] = mosaic_data['date']
  except: 
    fail('The Mosaic json (' + str(filename) + ') does not include a date of creation')

  # Check that the resource json reference matches the selected reference
  try:
    mosaic_info['reference'] = mosaic_data['reference']
  except:
    fail('The Mosaic json does not include a reference genome')
  if reference != mosaic_info['reference']:
    fail('The selected reference (' + str(reference) + ') does not match the Mosaic json reference (' + str(mosaic_info['reference']) + ')')

  # Store information on annotationns to add or remove from the project, and defaults to set
  if 'import_annotations' in mosaic_data:
    mosaic_info['import_annotations'] = mosaic_data['import_annotations'] if 'import_annotations' in mosaic_data else []
  if 'remove' in mosaic_data:
    mosaic_info['remove'] = mosaic_data['remove'] if 'remove' in mosaic_data else []
  if 'default_annotations' in mosaic_data:
    mosaic_info['defaultAnnotations'] = mosaic_data['default_annotations']

  # Store the name of the json describing the variant filters to be applied
  try:
    mosaic_info['variantFilters'] = mosaic_data['variant_filters']
  except:
    fail('The Mosaic json does not include a json describing the variant filters to apply')

  try:
    mosaic_info['exomiser_variant_filters'] = mosaic_data['exomiser_variant_filters']
  except:
    pass

  try:
    mosaic_info['sv_variant_filters'] = mosaic_data['sv_variant_filters']
  except:
    pass

  # Loop over all the specified resources, and get all information required to pull the annotations
  # into Mosaic
  try:
    resources = mosaic_data['resources']
  except:
    fail('The Mosaic json (' + str(filename) + ') does not include a list of resource')
  for resource in resources:
    mosaic_info['resources'][resource] = {}

    # The annotation must be identified as public or private. The other required information is dependent
    # on this
    try:
      mosaic_info['resources'][resource]['annotation_type'] = resources[resource]['annotation_type']
    except:
      fail('The Mosaic json file does not include the "annotation_type" field for resource: ' + str(resource))

    # The resource "class" dictates how the annotated vcf will be processed. If this is not present, it will be
    # processed in a standard way.
    mosaic_info['resources'][resource]['class'] = resources[resource]['class'] if 'class' in resources[resource] else False

    # Some annotations need to be extracted from an INFO field that is not the name provided in the annotation. For example, the HGVS
    # codes need to be extracted from the CSQ field. The "info_field" provides information on where the info should be extracted from
    mosaic_info['resources'][resource]['info_field'] = resources[resource]['info_field'] if 'info_field' in resources[resource] else False

    # Check if the "delimiter" value is set. This defines how to break up compound annotations
    mosaic_info['resources'][resource]['delimiter'] = resources[resource]['delimiter'] if 'delimiter' in resources[resource] else False

    # Collect the information required only for public annotations
    if mosaic_info['resources'][resource]['annotation_type'] == 'public':

      # Public annotations require the project id of the project that "owns" the annotation. Annotation values
      # need to be uploaded to the owner project
      try:
        project_id = resources[resource]['project_id']
      except:
        fail('The Mosaic json file does not include a "project_id" for public resource: ' + str(resource))
      mosaic_info['resources'][resource]['project_id'] = project_id

    # Fail if the supplied annotation_type is neither public or private
    elif mosaic_info['resources'][resource]['annotation_type'] != 'private':
      fail('Annotation_type for ' + str(resource) + ' must be "public" or "private"')

    # Loop over the annotation information
    mosaic_info['resources'][resource]['annotations'] = {}
    for annotation in resources[resource]['annotations']:
      mosaic_info['resources'][resource]['annotations'][annotation] = {}

      # Check that all required information is supplied and store it
      try:
        mosaic_info['resources'][resource]['annotations'][annotation]['uid'] = resources[resource]['annotations'][annotation]['uid']
      except:
        fail('The Mosaic json does not contain the "uid" field for annotation "' + str(annotation) + '" for resource "' + str(resource) + '"')
      try:
        mosaic_info['resources'][resource]['annotations'][annotation]['type'] = resources[resource]['annotations'][annotation]['type']
      except:
        fail('The Mosaic json does not contain the "type" field for annotation "' + str(annotation) + '" for resource "' + str(resource) + '"')

      # Check if the "is_sample" flag is set. This means there is an annotation for each project sample
      mosaic_info['resources'][resource]['annotations'][annotation]['isSample'] = resources[resource]['annotations'][annotation]['is_sample'] if 'is_sample' in resources[resource]['annotations'][annotation] else False

      # Check if the "position" value is set. This defines the position in a compound annotation that the desired annotation can be found
      mosaic_info['resources'][resource]['annotations'][annotation]['position'] = resources[resource]['annotations'][annotation]['position'] if 'position' in resources[resource]['annotations'][annotation] else False

      # Check if the "positions" value is set. This defines the position in a compound annotation that the desired annotation can be found
      mosaic_info['resources'][resource]['annotations'][annotation]['positions'] = resources[resource]['annotations'][annotation]['positions'] if 'positions' in resources[resource]['annotations'][annotation] else False

      # Check if the "operation" value is set. This indicates that the annotation is constructed from other annotations using a specific operation
      mosaic_info['resources'][resource]['annotations'][annotation]['operation'] = resources[resource]['annotations'][annotation]['operation'] if 'operation' in resources[resource]['annotations'][annotation] else False

      # Check if the "fields" value is set. This is used to determine which fields are used to generate the value for a specific operation
      mosaic_info['resources'][resource]['annotations'][annotation]['fields'] = resources[resource]['annotations'][annotation]['fields'] if 'fields' in resources[resource]['annotations'][annotation] else False

      # Check if "category" is set. If so, this will be used when creating new private annotations
      mosaic_info['resources'][resource]['annotations'][annotation]['category'] = resources[resource]['annotations'][annotation]['category'] if 'category' in resources[resource]['annotations'][annotation] else False

      # Check if "display_type" is set. If so, this will be used when creating new private annotations
      mosaic_info['resources'][resource]['annotations'][annotation]['display_type'] = resources[resource]['annotations'][annotation]['display_type'] if 'display_type' in resources[resource]['annotations'][annotation] else False

      # Check if "severity" is set. If so, the severity will be used when creating new private annotations
      mosaic_info['resources'][resource]['annotations'][annotation]['severity'] = resources[resource]['annotations'][annotation]['severity'] if 'severity' in resources[resource]['annotations'][annotation] else False

  # Return the mosaic information
  return mosaic_info

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = '')
  exit(1)
