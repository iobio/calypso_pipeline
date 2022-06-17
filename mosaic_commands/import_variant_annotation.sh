#!/bin/bash

TOKEN=$1
URL=$2
PROJECT_ID=$3
ATTRIBUTE_ID=$4

curl -S -s -X POST -H "Content-Type: application/json" -H "Authorization: Bearer $TOKEN" \
-d "{\"annotation_id\": $ATTRIBUTE_ID}" \
$URL"v1/projects/$PROJECT_ID/variants/annotations/import"
