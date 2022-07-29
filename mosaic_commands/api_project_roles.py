#!/usr/bin/python

# This contains API routes for project attributes (mirrors the API docs)

######
###### GET routes
######

# Get all roles attached to a project
def getProjectRoles(mosaicConfig, projectId, limit, page):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X GET -H "Authorization: Bearer ' + str(token) + '" '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/roles?limit=' + str(limit) + '&page=' + str(page)

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
