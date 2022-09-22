#!/usr/bin/python

from os.path import exists

# Parse the config file and get the Mosaic token and url
def parseConfig(args, required):
  token     = False
  url       = False
  projectId = False

  # Check the config file exists, if it was defined
  if args.config:
    if not exists(args.config): fail("Config file '" + str(args.config) + "' does not exist")

    # Parse the config file and store the token and url
    try: configFile = open(args.config, "r")
    except: fail("Failed to open config file '" + str(args.config) + "'")
    for line in configFile.readlines():
      fields = line.rstrip().split("=")
      if fields[0].startswith("MOSAIC_TOKEN"): token = fields[1]
      elif fields[0].startswith("MOSAIC_URL"): url = fields[1]
      elif fields[0].startswith("MOSAIC_ATTRIBUTES_PROJECT_ID"): projectId = fields[1]

  # Explicitly set attributes will overwrite the config file
  try: 
    if args.token: token = args.token
  except: token = False
  try:
    if args.url: url = args.url
  except: url = False
  try:
    if args.attributes_project: projectId = args.attributes_project
  except: projectId = False

  # Check that all required values are set
  if required["token"] and not token:
    print("An access token is required. You can either supply a token with '--token (-t)' or")
    print("supply a config file '--config (-c)' which includes the line:")
    print("  MOSAIC_TOKEN = <TOKEN>")
    exit(1)
  if required["url"] and not url:
    print("The api url is required. You can either supply the url with '--url (-u)' or")
    print("supply a config file '--config (-c)' which includes the line:")
    print("  MOSAIC_URL = <URL>")
    exit(1)
  if required["attributesProjectId"] and not projectId:
    print("The project id for the attributes project is required. You can either supply the id with '--attributesProject (-a)' or")
    print("supply a config file '--config (-c)' which includes the line:")
    print("  MOSAIC_ATTRIBUTES_PROJECT_ID = <ID>")
    exit(1)

  # Ensure the url terminates with a '/'
  if not url.endswith("/"): url += "/"

  # Return the values
  return {"token": token, "url": url, "attributesProjectId": projectId}

# If problems are found with the templates, fail
def fail(text):
  print(text)
  exit(1)
