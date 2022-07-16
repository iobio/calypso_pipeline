#!/usr/bin/python

# This contains API routes for samples (mirrors the API docs)

######
###### GET routes
######

# Get all samples in a project
def getSamples(mosaicConfig, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X GET -H "Authorization: Bearer ' + str(token) + '" '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/samples'

  return command

######
###### POST routes
######

######
###### PUT routes
######

######
###### DELETE routes
######
