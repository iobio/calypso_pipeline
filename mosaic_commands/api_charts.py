#!/usr/bin/python

# This contains API routes for charts (mirrors the API docs)

######
###### GET routes
######

# Get all project charts
def getProjectCharts(mosaicConfig, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X GET -H "Authorization: Bearer ' + str(token) + '" '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/charts'

  return command

######
###### POST routes
######

# Post a chart with background data
def postProjectBackgroundChart(mosaicConfig, name, chartType, attributeId, backgroundId, ylabel, colorId, compareId, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X POST -H "Content-Type: application/json" -H "Authorization: Bearer ' + str(token) + '" '
  command += '-d \'{"name": "' + str(name) + '", "chart_type": "' + str(chartType) + '", "attribute_id": "' + str(attributeId)
  command += '", "project_background_id": "' + str(backgroundId) + '", "saved_chart_data": {"x_axis": "attribute", "y_label": "'
  command += str(ylabel) + '", "color_by": "attribute", "color_by_attribute_id": "' + str(colorId) + '", "compare_to_attribute_id": "'
  command += str(compareId) + '"}}\' ' + str(url) + 'api/v1/projects/' + str(projectId) + '/charts'

  return command

######
###### PUT routes
######

######
###### DELETE routes
######

# Delete chart
def deleteSavedChart(mosaicConfig, projectId, chartId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X DELETE -H "Content-Type: application/json" -H "Authorization: Bearer ' + str(token) + '" '
  command += str(url) + 'api/v1/projects/' + str(projectId) + '/charts/' + str(chartId)

  return command
