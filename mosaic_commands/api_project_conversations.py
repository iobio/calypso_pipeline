#!/usr/bin/python

# This contains API routes for conversations (mirrors the API docs)

######
###### GET routes
######

# Get all project conversations
def getCoversations(mosaicConfig, limit, page, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X GET -H "Authorization: Bearer ' + str(token) + '" '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/conversations?limit=' + str(limit) + '&page=' + str(page)

  return command

######
###### POST routes
######

# Create a new project conversation
def postCoversation(mosaicConfig, title, description, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X POST -H "Content-Type: application/json" -H "Authorization: Bearer ' + str(token) + '" '
  command += '-d \'{"title": "' + str(title) + '", "description": "' + str(description) + '"}\' '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/conversations'

  return command

######
###### PUT routes
######

######
###### DELETE routes
######
