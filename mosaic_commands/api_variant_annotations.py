#!/usr/bin/python

# This contains API routes for variant annotations (mirrors the API docs)

######
###### GET routes
######

# Get the variant annotations in the project
def getVariantAnnotations(mosaicConfig, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X GET -H "Authorization: Bearer ' + str(token) + '" ' + str(url) + 'api/v1/projects/'
  command += str(projectId) + '/variants/annotations'

  return command

# Get variant annotations available to import
def getVariantAnnotationsImport(mosaicConfig, projectId, limit, page):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X GET -H "Authorization: Bearer ' + str(token) + '" ' + str(url) + 'api/v1/projects/'
  command += str(projectId) + '/variants/annotations/import?limit=' + str(limit) + '&page=' + str(page)

  return command

######
###### POST routes
######

# Create a variant annotation
def postCreateVariantAnnotations(mosaicConfig, name, valueType, privacy, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X POST -H "Content-Type: application/json" -H "Authorization: Bearer ' + str(token)
  command += '" -d \'{"name": "' + str(name) + '", "value_type": "' + str(valueType) + '", "privacy_level": "' + str(privacy) + '"}\' '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/variants/annotations'

  return command

# Import variant annotations to the project
def postImportVariantAnnotations(mosaicConfig, annotationId, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X POST -H "Content-Type: application/json" -H "Authorization: Bearer ' + str(token)
  command += '" -d \'{"annotation_id": ' + str(annotationId) + '}\' '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/variants/annotations/import'

  return command

######
###### PUT routes
######

######
###### DELETE routes
######

# Delete a variant annotation from a project
def deleteVariantAnnotation(mosaicConfig, projectId, annotationId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X DELETE -H "Authorization: Bearer ' + str(token) + '" '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/variants/annotations/' + str(annotationId)

  return command
