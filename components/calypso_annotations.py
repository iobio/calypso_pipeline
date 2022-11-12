#!/usr/bin/python

from __future__ import print_function

import os

#
def getAnnotationProjects(mosaicInfo, projectId):
  for resource in mosaicInfo['resources']:
    privacy       = mosaicInfo['resources'][resource]['annotation_type']
    annoProjectId = mosaicInfo['resources'][resource]['project_id'] if privacy == 'public' else projectId
    print(resource, privacy, annoProjectId)
  exit(0)

# Extract the annotations to be uploaded into Mosaic from the filtered vcf using bcftools
#def extractAnnotations(mosaicInfo, filteredVcf):

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
