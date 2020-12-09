#!/bin/bash


# Assign taxonomy for OTUs
python scripts/assign_taxonomy_for_consensus_seqs.py \
    --rdp-table files/assign_taxonomy_OTUs/taxonomy_RDP.txt \
    --confidence-threshold 50 \
    --output-basepath ../output/processed_data

# Filter the OTUs using consistency filtering
python ../step_2_filtering.py \
    --dataset ../processed_data/gibson_healthy_agg_taxa.pkl \
    --outfile ../processed_data/gibson_healthy_agg_taxa_filtered.pkl \
    --dtype rel \
    --threshold 0.0001 \
    --min-num-consecutive 7 \
    --min-num-subjects 2 \
    --colonization-time 5
python ../step_2_filtering.py \
    --dataset ../processed_data/gibson_uc_agg_taxa.pkl \
    --outfile ../processed_data/gibson_uc_agg_taxa_filtered.pkl \
    --dtype rel \
    --threshold 0.0001 \
    --min-num-consecutive 7 \
    --min-num-subjects 2 \
    --colonization-time 5

# Learn negative binomial dispersion parameters
python scripts/filter_replicates_like_other_dataset.py \
    --replicate-dataset ../processed_data/gibson_replicates_agg_taxa.pkl \
    --like-other ../processed_data/gibson_healthy_agg_taxa_filtered.pkl \
    --output-basepath ../processed_data/gibson_replicates_agg_taxa_filtered.pkl
