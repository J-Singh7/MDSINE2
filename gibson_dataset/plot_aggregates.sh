#!/bin/bash

# Plot the OTU aggregates
python scripts/plot_otus.py \
    --study ../processed_data/gibson_healthy_agg_taxa.pkl \
    --output-basepath ../output/aggregate_plots/
python scripts/plot_otus.py \
    --study ../processed_data/gibson_uc_agg_taxa.pkl \
    --output-basepath ../output/aggregate_plots/

