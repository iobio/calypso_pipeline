from os.path import exists

import os

# Build the toml file defining all the annotations that vcfanno should use
def buildToml(workingDir, resourceInfo):

  # Create a toml file
  tomlFilename  = 'calypso_annotations.toml'
  toml          = str(workingDir) + str(tomlFilename)
  try: tomlFile = open(toml, 'w')
  except: fail('There was a problem opening a file (calypso_annotations.toml) to write to')

  # Get the path to the data
  datapath = resourceInfo['path']

  # Add each required resource to the toml file
  for resource in resourceInfo['resources']:
    if resourceInfo['resources'][resource]['toml']:
      print('[[annotation]]', file = tomlFile)
      print('file="', str(datapath), str(resourceInfo['resources'][resource]['file']), '"', sep = '', file = tomlFile)

      # Some of the annotations include "fields", "columns", and "names". Include these in the toml if they exist
      if 'fields' in resourceInfo['resources'][resource]:
        text, noValues = tomlInfo(resourceInfo, resource, 'fields')
        print(text, file = tomlFile)
      if 'columns' in resourceInfo['resources'][resource]:
        text, noValues = tomlInfo(resourceInfo, resource, 'columns')
        print(text, file = tomlFile)
      if 'names' in resourceInfo['resources'][resource]:
        text, noValues = tomlInfo(resourceInfo, resource, 'names')
        print(text, file = tomlFile)

      # Write out the ops
      ops   = 'ops=["'
      noOps = len(resourceInfo['resources'][resource]['ops'])
      for i in range(0, noOps - 1): ops += str(resourceInfo['resources'][resource]['ops'][i]) + '", "'
      ops += str(resourceInfo['resources'][resource]['ops'][-1]) + '"]'
      print(ops, file = tomlFile)
      print(file = tomlFile)

      # If there is post annotation information to include in the toml, include it
      if 'post_annotation' in resourceInfo['resources'][resource]:
        post = resourceInfo['resources'][resource]['post_annotation']
        print('[[postannotation]]', file = tomlFile)
        if 'name' in post: print('name="', post['name'], '"', sep = '', file = tomlFile)
        if 'fields' in post:
          fieldString = 'fields=["'
          for i in range(0, len(post['fields']) - 1): fieldString += str(post['fields'][i]) + '", "'
          fieldString += str(post['fields'][-1]) + '"]'
          print(fieldString, file = tomlFile)
        if 'op' in post: print('op="', post['op'], '"', sep = '', file = tomlFile)
        if 'type' in post: print('type="', post['type'], '"', sep = '', file = tomlFile)
        print(file = tomlFile)

  # Close the toml file
  tomlFile.close()

  return tomlFilename

# Include information in the toml file
def tomlInfo(resourceInfo, resource, infoType):
  noValues = len(resourceInfo['resources'][resource][infoType])

  # Define the fields to use
  text = infoType + '=['
  for index, value in enumerate(resourceInfo['resources'][resource][infoType]):
    if infoType == 'columns': text += str(value)
    else: text += '"' + str(value) + '"'
    if (index + 1) < noValues: text += ', '
  text += ']'

  return text, noValues

# Generate a lua file with required functions
def generateLuaFile(workingDir):

  # Create script file for running exomiser
  scriptBase = 'calypso_vcfanno_lua.lua'
  scriptName = str(workingDir) + str(scriptBase)
  script     = open(scriptName, 'w')
  print('function hemi(nonpar, xy)', file = script)
  print('  if (nonpar == true)', file = script)
  print('  then', file = script)
  print('    return string.format("%d", xy)', file = script)
  print('  end', file = script)
  print('end', file = script)

  # Close the exomiser script file
  script.close()

  # Return the name of the lua file
  return scriptBase

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
