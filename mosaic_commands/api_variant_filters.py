#!/usr/bin/python

import json

# This contains API routes for variant filters (mirrors the API docs)

######
###### GET routes
######

# Get the variant filters in the project
def getVariantFilters(mosaicConfig, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X GET -H "Authorization: Bearer ' + str(token) + '" ' + str(url) + 'api/v1/projects/'
  command += str(projectId) + '/variants/filters'

  return command

######
###### POST routes
######

# Post a new variant filter
def postVariantFilter(mosaicConfig, name, annotationFilters, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X POST -H "Content-Type: application/json" -H "Authorization: Bearer ' + str(token) + '" '
  command += '-d \'{"name": "' + str(name) + '", "filter": ' + json.dumps(annotationFilters) + '}\' '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/variants/filters'

  return command

######
###### PUT routes
######

# Update a variant filter
def putVariantFilter(mosaicConfig, name, annotationFilters, projectId, filterId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X PUT -H "Content-Type: application/json" -H "Authorization: Bearer ' + str(token) + '" '
  command += '-d \'{"name": "' + str(name) + '", "filter": ' + json.dumps(annotationFilters) + '}\' '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/variants/filters/' + str(filterId)

  return command

######
###### DELETE routes
######

# Delete a variant filter
def deleteVariantFilter(mosaicConfig, projectId, filterId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X DELETE -H "Authorization: Bearer ' + str(token) + '" '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/variants/filters/' + str(filterId)

  return command
