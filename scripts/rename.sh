#!/bin/bash

while IFS=, read -r SRR NEWNAME; do
    echo "Renaming $SRR to $NEWNAME..."

    mv raw_data/"${SRR}_1.fastq.gz" raw_data/"${NEWNAME}-R1.fastq.gz"
    mv raw_data/"${SRR}_2.fastq.gz" raw_data/"${NEWNAME}-R2.fastq.gz"
done < config/mapping.csv