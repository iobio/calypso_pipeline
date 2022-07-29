#!/usr/bin/python

# This contains API routes for dashboards (mirrors the API docs)

######
###### GET routes
######

# Get dashboard information
def getDashboard(mosaicConfig, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X GET -H "Authorization: Bearer ' + str(token) + '" '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/dashboard'

  return command

######
###### POST routes
######

# Pin a chart to the dashboard
def postPinChart(mosaicConfig, chartId, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X POST -H "Content-Type: application/json" -H "Authorization: Bearer ' + str(token) + '" '
  command += '-d \'{"type": "chart", "chart_id": "' + str(chartId) + '", "is_active": "true"}\' '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/dashboard'

  return command

# Pin an attribute to the dashboard
def postPinAttribute(mosaicConfig, attributeId, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X POST -H "Content-Type: application/json" -H "Authorization: Bearer ' + str(token) + '" '
  command += '-d \'{"type": "project_attribute", "attribute_id": "' + str(attributeId) + '", "is_active": "true"}\' '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/dashboard'

  return command

# Pin a conversation to the dashboard
def postPinConversation(mosaicConfig, conversationId, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X POST -H "Content-Type: application/json" -H "Authorization: Bearer ' + str(token) + '" '
  command += '-d \'{"type": "conversation", "project_conversation_id": "' + str(conversationId) + '", "is_active": "true"}\' '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/dashboard'

  return command

######
###### PUT routes
######

######
###### DELETE routes
######
