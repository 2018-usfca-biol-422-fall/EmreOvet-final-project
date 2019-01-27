#!/bin/bash

# the script to download SRA files from NCBI

# Emre Ovet
# November 28, 2018
# eovet@dons.usfca.edu

# Make sure to download SRA Toolkit online before running this code
for SRA_number in $(cut -f 6 data/metadata/SraRunTable.txt | tail -n +2)
do
  /Users/emreovet/sratoolkit.2.9.2-mac64/bin/fastq-dump -v "$SRA_number" -O data/raw_data
done
