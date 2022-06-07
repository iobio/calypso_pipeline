#!/bin/bash

# TYPE="allele | position"
# NAME="Name of variant set"

TOKEN=$1
URL=$2
PROJECT_ID=$3
FILENAME=$4
TYPE=$5
NAME=$6

curl -S -s -i -X POST -H "Content-Type: multipart/form-data" -H "Authorization: Bearer $TOKEN" \
-F "file=@$FILENAME" \
-F "type=$TYPE" \
-F "name=$NAME" \
$URL"v1/projects/$PROJECT_ID/variants/upload"
