#!/bin/bash

# Agglomerate ASVs into OTUs
python scripts/preprocess.py \
    --hamming-distance 2 \
    --rename-prefix OTU \
    --sequences files/preprocessing/gibson_16S_rRNA_v4_ASV_seqs_aligned_filtered.fa \
    --output-basepath ../output/processed_data \
    --remove-timepoints 0 0.5