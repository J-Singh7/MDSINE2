#!/bin/bash

# Plot the phylogenetic subtrees for each OTU
python scripts/plot_phylogenetic_subtrees.py \
    --study ../processed_data/gibson_healthy_agg_taxa.pkl \
    --output-basepath ../output/subtrees/ \
    --tree files/phylogenetic_placement_OTUs/phylogenetic_tree_full_taxid.nhx \
    --seq-info files/subtrees/RDP-11-5_BA_TS_info.tsv \
    --sep \t \
    --family-radius-factor 1.5
