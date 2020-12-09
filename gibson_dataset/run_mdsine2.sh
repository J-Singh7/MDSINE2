#!/bin/bash

NEGBIN="../output/negbin/replicates/mcmc.pkl"
BURNIN="5000"
N_SAMPLES="15000"
CHECKPOINT="100"
MULTIPROCESSING="1"
HEALTHY_DSET="../processed_data/gibson_healthy_agg_taxa_filtered.pkl"
UC_DSET="../processed_data/gibson_uc_agg_taxa_filtered.pkl"
INTERACTION_IND_PRIOR="strong-sparse"
PERTURBATION_IND_PRIOR="strong-sparse"

# Healthy cohort
# --------------
python ../step_5_infer_mdsine2.py \
    --input $HEALTHY_DSET \
    --negbin $NEGBIN \
    --seed 0 \
    --burnin $BURNIN \
    --n-samples $N_SAMPLES \
    --checkpoint $CHECKPOINT \
    --multiprocessing 1 \
    --rename-study healthy-seed0 \
    --basepath ../output/mdsine2 \
    --interaction-ind-prior $INTERACTION_IND_PRIOR \
    --perturbation-ind-prior $PERTURBATION_IND_PRIOR
python ../step_6_visualize_mdsine2.py \
    --chain  ../output/mdsine2/healthy-seed0/mcmc.pkl
    --output-basepath ../output/mdsine2/healthy-seed0/posterior

python ../step_5_infer_mdsine2.py \
    --input $HEALTHY_DSET \
    --negbin $NEGBIN \
    --seed 1 \
    --burnin $BURNIN \
    --n-samples $N_SAMPLES \
    --checkpoint $CHECKPOINT \
    --multiprocessing 1 \
    --rename-study healthy-seed1 \
    --basepath ../output/mdsine2 \
    --interaction-ind-prior $INTERACTION_IND_PRIOR \
    --perturbation-ind-prior $PERTURBATION_IND_PRIOR
python ../step_6_visualize_mdsine2.py \
    --chain  ../output/mdsine2/healthy-seed1/mcmc.pkl
    --output-basepath ../output/mdsine2/healthy-seed1/posterior

# UC cohort
# ---------
python ../step_5_infer_mdsine2.py \
    --input $UC_DSET \
    --negbin $NEGBIN \
    --seed 0 \
    --burnin $BURNIN \
    --n-samples $N_SAMPLES \
    --checkpoint $CHECKPOINT \
    --multiprocessing 1 \
    --rename-study uc-seed0 \
    --basepath ../output/mdsine2 \
    --interaction-ind-prior $INTERACTION_IND_PRIOR \
    --perturbation-ind-prior $PERTURBATION_IND_PRIOR
python ../step_6_visualize_mdsine2.py \
    --chain  ../output/mdsine2/uc-seed0/mcmc.pkl
    --output-basepath ../output/mdsine2/uc-seed0/posterior

python ../step_5_infer_mdsine2.py \
    --input $UC_DSET \
    --negbin $NEGBIN \
    --seed 1 \
    --burnin $BURNIN \
    --n-samples $N_SAMPLES \
    --checkpoint $CHECKPOINT \
    --multiprocessing 1 \
    --rename-study uc-seed0 \
    --basepath ../output/mdsine2 \
    --interaction-ind-prior $INTERACTION_IND_PRIOR \
    --perturbation-ind-prior $PERTURBATION_IND_PRIOR
python ../step_6_visualize_mdsine2.py \
    --chain  ../output/mdsine2/uc-seed1/mcmc.pkl
    --output-basepath ../output/mdsine2/uc-seed1/posterior
