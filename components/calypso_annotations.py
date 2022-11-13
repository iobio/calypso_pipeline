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
def createAnnotationTsv(mosaicInfo, resource, scriptDir, configFile, vcf, privateAnnotations, bashFile):
  tags = ''
  uids = ''

  # Public resources will have a uid associated with the resource and so this can be extracted
  for annotation in mosaicInfo['resources'][resource]['annotations']:
    tags += ',' + annotation
    uids += ',' + mosaicInfo['resources'][resource]['annotations'][annotation]['uid']
  tags = tags.strip(',')
  uids = uids.strip(',')

  print(file = bashFile)
  print('  # Resource: ', resource, sep = '', file = bashFile)
  print('  echo -n "Creating tsv file for resource ', resource, '..."', sep = '', file = bashFile)
  print('  python ', scriptDir, '/generate_annotation_tsv.py \\', sep = '', file = bashFile)
  print('    -c "', configFile, '" \\', sep = '', file = bashFile)
  print('    -r "', resource, '" \\', sep = '', file = bashFile)
  print('    -g "', tags, '" \\', sep = '', file = bashFile)
  print('    -d "', uids, '" \\', sep = '', file = bashFile)
  print('    -i "', vcf, '"', sep = '', file = bashFile)
  print('  echo "complete"', file = bashFile)

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
