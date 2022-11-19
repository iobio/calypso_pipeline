#!/usr/bin/python

from __future__ import print_function
from datetime import date
from os.path import exists
from sys import path

import os
import argparse
import json
import math
import glob
import importlib

# Compare the current and previous resources and identify any differences
def compare(previousResourceInfo, resourceInfo):
  addedResources   = []
  removedResources = list(previousResourceInfo['resources'].keys())
  updatedResources = {}

  # Loop over the resources in the current resource json and compare to the previous json
  for resource in resourceInfo['resources']:
    if resource not in previousResourceInfo['resources']: addedResources.append(resource)
    else:
      if str(resourceInfo['resources'][resource]['version']) != str(previousResourceInfo['resources'][resource]['version']):
        updatedResources[resource] = {'previousVersion': previousResourceInfo['resources'][resource]['version'], 'currentVersion': resourceInfo['resources'][resource]['version']}
      removedResources.remove(resource)

  # Write out a summary of changes
  if len(removedResources) > 0:
    print()
    print('The following resources have been removed from Calypso so will not be updated:')
    for resource in removedResources: print('  ', resource)
  if len(addedResources) > 0:
    print()
    print('The following resources have been added to Calypso, so this will be the first inclusion of these annotations:')
    for resource in addedResources: print('  ', resource)
  if len(updatedResources) > 0:
    print()
    print('The following resources have been updated:')
    for resource in updatedResources: print('  ', resource, ' was updated from ', updatedResources[resource]['previousVersion'], ' to ', updatedResources[resource]['currentVersion'], sep = '')

  # Return the updated resources
  return updatedResources

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
