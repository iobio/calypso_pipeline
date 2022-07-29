#!/usr/bin/python

# This contains API routes for project attributes (mirrors the API docs)

######
###### GET routes
######

######
###### POST routes
######

# Create a project
def postProject(mosaicConfig, name, reference):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X POST -H "Content-Type: application/json" -H "Authorization: Bearer ' + str(token) + '" '
  command += '-d \'{"name": "' + str(name) + '", "reference": "' + str(reference) + '"}\' ' + str(url) + 'api/v1/projects'

  return command

######
###### PUT routes
######

######
###### DELETE routes
######

