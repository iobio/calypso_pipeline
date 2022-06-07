#!/bin/bash

TOKEN=$1
URL=$2
PROJECT_ID=$3
FILENAME=$4

curl -s -S -i -X POST -H "Content-Type: multipart/form-data" -H "Authorization: Bearer $TOKEN" \
-F "file=@$FILENAME" \
$URL"v1/projects/$PROJECT_ID/variants/annotations/upload"
