#!/usr/bin/python

# This contains API routes for project backgrounds (mirrors the API docs)

######
###### GET routes
######

######
###### POST routes
######

# Post background data to a project
def postBackgrounds(mosaicConfig, filename, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X POST -H "Content-Type: application/json" -H "Authorization: Bearer ' + str(token) + '" '
  command += '-d "@' + str(filename) + '" ' + str(url) + 'api/v1/projects/' + str(projectId) + '/backgrounds'

  return command

# Post a background chart to a project
def postBackgroundChart(mosaicConfig, name, chartType, attributeId, backgroundId, yLabel, colour, compare, projectId):
  token = mosaicConfig["token"]
  url   = mosaicConfig["url"]

  command  = 'curl -S -s -X POST -H "Content-Type: application/json" -H "Authorization: Bearer ' + str(token) + '" '
  command += '-d \'{"name": "' + str(name) + '", "chart_type": "' + str(chartType) + '", "attribute_id": "' + str(attributeId)
  command += '", "project_background_id": "' + str(backgroundId) +'", "saved_chart_data": {"x_axis": "attribute", "y_label": "'
  command += str(yLabel) + '", "color_by": "attribute", "color_by_attribute_id": ' + str(colour) + ', "compare_to_attribute_id": '
  command += str(compare) + '}}\' ' + str(url) + 'api/v1/projects/' + str(projectId) + '/charts'

  return command

######
###### PUT routes
######

######
###### DELETE routes
######
