#!/bin/bash

NEGBIN="../output/negbin/replicates/mcmc.pkl"
SEED="0"
BURNIN="5000"
N_SAMPLES="15000"
CHECKPOINT="100"
MULTIPROCESSING="0"
DSET_BASEPATH="../processed_data/cv"
CV_BASEPATH="../output/mdsine2/cv"
HEALTHY_DATASET="../processed_data/gibson_healthy_agg_taxa_filtered.pkl"
UC_DATASET="../processed_data/gibson_uc_agg_taxa_filtered.pkl"
INTERACTION_IND_PRIOR="weak-agnostic"
PERTURBATION_IND_PRIOR="weak-agnostic"

# Run healthy for each subject
# ----------------------------
# Regular inference
python ../run_cross_validation.py \
    --dataset $HEALTHY_DATASET \
    --cv-basepath $CV_BASEPATH \
    --dset-basepath $DSET_BASEPATH \
    --negbin $NEGBIN \
    --seed $SEED \
    --burnin $BURNIN \
    --n-samples $N_SAMPLES \
    --checkpoint $CHECKPOINT \
    --multiprocessing $MULTIPROCESSING \
    --leave-out-subject 2 \
    --interaction-ind-prior $INTERACTION_IND_PRIOR \
    --perturbation-ind-prior $PERTURBATION_IND_PRIOR
CMD="python ../step_6_visualize_mdsine2.py \
    --chain ${CV_BASEPATH}/healthy-cv2/mcmc.pkl
    --output-basepath ${CV_BASEPATH}/healthy-cv2/posterior"
eval $CMD

python ../run_cross_validation.py \
    --dataset $HEALTHY_DATASET \
    --cv-basepath $CV_BASEPATH \
    --dset-basepath $DSET_BASEPATH \
    --negbin $NEGBIN \
    --seed $SEED \
    --burnin $BURNIN \
    --n-samples $N_SAMPLES \
    --checkpoint $CHECKPOINT \
    --multiprocessing $MULTIPROCESSING \
    --leave-out-subject 3 \
    --interaction-ind-prior $INTERACTION_IND_PRIOR \
    --perturbation-ind-prior $PERTURBATION_IND_PRIOR
CMD="python ../step_6_visualize_mdsine2.py \
    --chain ${CV_BASEPATH}/healthy-cv3/mcmc.pkl
    --output-basepath ${CV_BASEPATH}/healthy-cv3/posterior"
eval $CMD

python ../run_cross_validation.py \
    --dataset $HEALTHY_DATASET \
    --cv-basepath $CV_BASEPATH \
    --dset-basepath $DSET_BASEPATH \
    --negbin $NEGBIN \
    --seed $SEED \
    --burnin $BURNIN \
    --n-samples $N_SAMPLES \
    --checkpoint $CHECKPOINT \
    --multiprocessing $MULTIPROCESSING \
    --leave-out-subject 4 \
    --interaction-ind-prior $INTERACTION_IND_PRIOR \
    --perturbation-ind-prior $PERTURBATION_IND_PRIOR
CMD="python ../step_6_visualize_mdsine2.py \
    --chain ${CV_BASEPATH}/healthy-cv4/mcmc.pkl
    --output-basepath ${CV_BASEPATH}/healthy-cv4/posterior"
eval $CMD

python ../run_cross_validation.py \
    --dataset $HEALTHY_DATASET \
    --cv-basepath $CV_BASEPATH \
    --dset-basepath $DSET_BASEPATH \
    --negbin $NEGBIN \
    --seed $SEED \
    --burnin $BURNIN \
    --n-samples $N_SAMPLES \
    --checkpoint $CHECKPOINT \
    --multiprocessing $MULTIPROCESSING \
    --leave-out-subject 5 \
    --interaction-ind-prior $INTERACTION_IND_PRIOR \
    --perturbation-ind-prior $PERTURBATION_IND_PRIOR
CMD="python ../step_6_visualize_mdsine2.py \
    --chain ${CV_BASEPATH}/healthy-cv5/mcmc.pkl
    --output-basepath ${CV_BASEPATH}/healthy-cv5/posterior"
eval $CMD

# Run uc for each subject
# -----------------------
python ../run_cross_validation.py \
    --dataset $uc_DATASET \
    --cv-basepath $CV_BASEPATH \
    --dset-basepath $DSET_BASEPATH \
    --negbin $NEGBIN \
    --seed $SEED \
    --burnin $BURNIN \
    --n-samples $N_SAMPLES \
    --checkpoint $CHECKPOINT \
    --multiprocessing $MULTIPROCESSING \
    --leave-out-subject 6 \
    --interaction-ind-prior $INTERACTION_IND_PRIOR \
    --perturbation-ind-prior $PERTURBATION_IND_PRIOR
CMD="python ../step_6_visualize_mdsine2.py \
    --chain ${CV_BASEPATH}/uc-cv6/mcmc.pkl
    --output-basepath ${CV_BASEPATH}/uc-cv6/posterior"
eval $CMD

python ../run_cross_validation.py \
    --dataset $uc_DATASET \
    --cv-basepath $CV_BASEPATH \
    --dset-basepath $DSET_BASEPATH \
    --negbin $NEGBIN \
    --seed $SEED \
    --burnin $BURNIN \
    --n-samples $N_SAMPLES \
    --checkpoint $CHECKPOINT \
    --multiprocessing $MULTIPROCESSING \
    --leave-out-subject 7 \
    --interaction-ind-prior $INTERACTION_IND_PRIOR \
    --perturbation-ind-prior $PERTURBATION_IND_PRIOR
CMD="python ../step_6_visualize_mdsine2.py \
    --chain ${CV_BASEPATH}/uc-cv7/mcmc.pkl
    --output-basepath ${CV_BASEPATH}/uc-cv7/posterior"
eval $CMD

python ../run_cross_validation.py \
    --dataset $uc_DATASET \
    --cv-basepath $CV_BASEPATH \
    --dset-basepath $DSET_BASEPATH \
    --negbin $NEGBIN \
    --seed $SEED \
    --burnin $BURNIN \
    --n-samples $N_SAMPLES \
    --checkpoint $CHECKPOINT \
    --multiprocessing $MULTIPROCESSING \
    --leave-out-subject 8 \
    --interaction-ind-prior $INTERACTION_IND_PRIOR \
    --perturbation-ind-prior $PERTURBATION_IND_PRIOR
CMD="python ../step_6_visualize_mdsine2.py \
    --chain ${CV_BASEPATH}/uc-cv8/mcmc.pkl
    --output-basepath ${CV_BASEPATH}/uc-cv8/posterior"
eval $CMD

python ../run_cross_validation.py \
    --dataset $uc_DATASET \
    --cv-basepath $CV_BASEPATH \
    --dset-basepath $DSET_BASEPATH \
    --negbin $NEGBIN \
    --seed $SEED \
    --burnin $BURNIN \
    --n-samples $N_SAMPLES \
    --checkpoint $CHECKPOINT \
    --multiprocessing $MULTIPROCESSING \
    --leave-out-subject 9 \
    --interaction-ind-prior $INTERACTION_IND_PRIOR \
    --perturbation-ind-prior $PERTURBATION_IND_PRIOR
CMD="python ../step_6_visualize_mdsine2.py \
    --chain ${CV_BASEPATH}/uc-cv9/mcmc.pkl
    --output-basepath ${CV_BASEPATH}/uc-cv9/posterior"
eval $CMD

python ../run_cross_validation.py \
    --dataset $uc_DATASET \
    --cv-basepath $CV_BASEPATH \
    --dset-basepath $DSET_BASEPATH \
    --negbin $NEGBIN \
    --seed $SEED \
    --burnin $BURNIN \
    --n-samples $N_SAMPLES \
    --checkpoint $CHECKPOINT \
    --multiprocessing $MULTIPROCESSING \
    --leave-out-subject 10 \
    --interaction-ind-prior $INTERACTION_IND_PRIOR \
    --perturbation-ind-prior $PERTURBATION_IND_PRIOR
CMD="python ../step_6_visualize_mdsine2.py \
    --chain ${CV_BASEPATH}/uc-cv10/mcmc.pkl
    --output-basepath ${CV_BASEPATH}/uc-cv10/posterior"
eval $CMD