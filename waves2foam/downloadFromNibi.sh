#!/bin/bash
# Script to download simulation results from Nibi to local workstation
# Usage: ./downloadFromNibi.sh

# Your Nibi username
USER=andyhe

# Path on Nibi where the case is
CASE_PATH=/scratch/andyhe/work/wavemakerTank

# Name for the local folder
CASE_NAME=wavemakerTank

# Download using rsync
rsync -av --no-g --no-p \
    ${USER}@nibi.alliancecan.ca:${CASE_PATH} \
    ./${CASE_NAME}

echo "Download complete!"
