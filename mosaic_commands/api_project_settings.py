#!/usr/bin/python

# This contains API routes for project settings (mirrors the API docs)

######
###### GET routes
######

# Get the project settings
def getProjectSettings(mosaicConfig, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command = 'curl -S -s -X GET -H "Authorization: Bearer ' + str(token) + '" ' + str(url) + 'api/v1/projects/' + str(projectId) + '/settings'

  return command

######
###### POST routes
######

######
###### PUT routes
######

# Set a variant annotation as a default
def putDefaultAnnotations(mosaicConfig, projectId, annotationIds):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X PUT -H "Content-Type: application/json" -H "Authorization: Bearer ' + str(token) + '" '
  command += '-d \'{"selected_variant_annotation_ids": ' + str(annotationIds) + '}\' '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/settings'

  return command

######
###### DELETE routes
######
