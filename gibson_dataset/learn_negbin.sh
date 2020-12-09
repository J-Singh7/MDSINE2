#!/bin/bash


python ../step_3_infer_negbin.py \
    --input ../processed_data/gibson_replicates_agg_taxa_filtered.pkl \
    --seed 0 \
    --burnin 2000 \
    --n-samples 6000 \
    --checkpoint 200 \
    --multiprocessing 0 \
    --basepath ../output/negbin/

python ../step_4_visualize_negbin.py \
    --chain ../output/negbin/replicates/mcmc.pkl \
    --output-basepath ../output/negbin/replicates/posterior

