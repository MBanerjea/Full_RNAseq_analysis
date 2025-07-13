#!/bin/bash

# Set number of threads
THREADS=4

# Directory to save downloads
mkdir -p raw_data

# Loop through each SRR ID
for SRR in $(cat config/SRR_Acc_List.txt); do
    echo "Downloading $SRR..."
    
    # Download paired-end reads using fasterq-dump
    prefetch ${SRR}
    echo "Splitting ${SRR}..."
    fasterq-dump --split-files -e $THREADS ${SRR} -O raw_data/

    # Gzip compress the FASTQ files
    echo "Compressing ${SRR}..."
    pigz -p ${THREADS} raw_data/${SRR}_1.fastq
    pigz -p ${THREADS} raw_data/${SRR}_2.fastq

    # Optional: rename to something human-readable
    # e.g., WT_1A_R1.fastq.gz — depending on your naming strategy

done