#!/usr/bin/python

# This contains API routes for attributes (mirrors the API docs)

######
###### GET routes
######

# Get all available attribute sets
def getProjectAttributeSets(mosaicConfig, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X GET -H "Authorization: Bearer ' + str(token) + '" '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/attributes/sets'

  return command

######
###### POST routes
######

# Create an attribute set. ValueType should be one of "project", "sample", or "gene"
def postProjectAttributeSet(mosaicConfig, name, description, isPublic, attributeIds, valueType, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X POST -H "Content-Type: application/json" -H "Authorization: Bearer ' + str(token) + '" '
  command += '-d \'{"name": "' + str(name) + '", "description": "' + str(description) + '", "is_public_to_project": "' + str(isPublic)
  command += '", "attribute_ids": "' + str(attributeIds) + '", "type": "' + str(valueType) + '"}\' '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/attributes/sets'

  return command

######
###### PUT routes
######

######
###### DELETE routes
######

# Delete an attribute set
def deleteProjectAttributeSet(mosaicConfig, projectId, setId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X DELETE -H "Content-Type: application/json" -H "Authorization: Bearer ' + str(token) + '" '
  command += url + 'api/v1/projects/' + str(projectId) + '/attributes/sets/' + str(setId)

  return command
