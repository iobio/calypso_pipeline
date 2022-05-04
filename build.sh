#!/bin/bash

set -eou pipefail

sudo docker build -t calypso_pipeline .
sudo apptainer build calypso_pipeline.sif docker-daemon://calypso_pipeline:latest
