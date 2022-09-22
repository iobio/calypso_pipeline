#!/usr/bin/python

import json

# This contains API routes for variant (mirrors the API docs)

######
###### GET routes
######

######
###### POST routes
######

# Upload variant via a file
def postUploadVariants(mosaicConfig, filename, fileType, name, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X POST -H "Content-Type: multipart/form-data" -H "Authorization: Bearer ' + str(token) + '" '
  command += '-F "file=@' + str(filename) + '" -F "type=' + str(fileType) + '" -F "name=' + str(name) + '" '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/variants/upload'

  return command

######
###### PUT routes
######

######
###### DELETE routes
######

