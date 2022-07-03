#!/usr/bin/python

# This contains API routes for project attributes (mirrors the API docs)

######
###### GET routes
######

# Get the project attributes for the defined project
def getProjectAttributes(mosaicConfig, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command = 'curl -S -s -X GET -H "Authorization: Bearer ' + str(token) + '" ' + str(url) + 'api/v1/projects/' + str(projectId) + '/attributes'

  return command

# Get all public project attributes
def getPublicProjectAttributes(mosaicConfig, limit, page):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X GET -H "Authorization: Bearer ' + str(token) + '" '
  command += str(url) + 'api/v1/projects/attributes?limit=' + str(limit) + '&page=' + str(page)

  return command

# Get the project attributes for the all projects the user has access to
def getUserProjectAttributes(mosaicConfig):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command = 'curl -S -s -X GET -H "Authorization: Bearer ' + str(token) + '" ' + str(url) + 'api/v1/user/projects/attributes'

  return command

######
###### POST routes
######

# Create a new project attribute
def postProjectAttribute(mosaicConfig, name, valueType, value, isPublic, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X POST -H "Content-Type: application/json" -H "Authorization: Bearer ' + str(token) + '" '
  command += '-d \'{"name": "' + str(name) + '", "value_type": "' + str(valueType) + '", "value": "' + str(value) + '", "is_public": "' + str(isPublic) + '"}\' '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/attributes'

  return command

# Import a project attribute
def postImportProjectAttribute(mosaicConfig, attributeId, value, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X POST -H "Content-Type: application/json" -H "Authorization: Bearer ' + str(token) + '" '
  command += '-d \'{"attribute_id": "' + str(attributeId) + '", "value": "' + str(value) + '"}\' '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/attributes/import'

  return command

######
###### PUT routes
######

# Update the value of a project attribute
def putProjectAttribute(mosaicConfig, value, projectId, attributeId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X PUT -H "Content-Type: application/json" -H "Authorization: Bearer ' + str(token) + '" '
  command += '-d \'{"value": "' + str(value) + '"}\' '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/attributes/' + str(attributeId)

  return command

######
###### DELETE routes
######
