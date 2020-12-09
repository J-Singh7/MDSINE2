#!/bin/bash

# ErisOne parameters
# ------------------
# Path to MDSINE2_Paper code
MDSINE2_PAPER_CODE_PATH="/data/cctm/darpa_perturbation_mouse_study/MDSINE2_Paper"
# Conda environment
ENVIRONMENT_NAME="mdsine2"
# Queues, memory, and numpy of cpus
QUEUE="short"
MEM="4000"
N_CPUS="1"
LSF_BASEPATH="lsf_files/tla"

# Forward simulation parameters
# -----------------------------
SIM_DT="0.01"
N_DAYS="8"
OUTPUT_BASEPATH="output/mdsine2/cv/forward_sims"

# Healthy cohort
# --------------
python scripts/run_forward_sim_for_fold.py \
    --chain output/mdsine2/cv/healthy-cv2/numpy_trace \
    --validation processed_data/cv/healthy-cv2-validate.pkl \
    --validation-curr-path ../../processed_data/cv/healthy-cv2-validate.pkl \
    --simulation-dt $SIM_DT \
    --n-days $N_DAYS \
    --output-basepath $OUTPUT_BASEPATH \
    --environment-name $ENVIRONMENT_NAME \
    --code-basepath $MDSINE2_PAPER_CODE_PATH \
    --queue $QUEUE \
    --memory $MEM \
    --n-cpus $N_CPUS \
    --lsf-basepath $LSF_BASEPATH

# python scripts/run_forward_sim_for_fold.py \
#     --chain output/mdsine2/cv/healthy-cv3/numpy_trace \
#     --validation processed_data/cv/healthy-cv3-validate.pkl \
#     --validation-curr-path ../../processed_data/cv/healthy-cv3-validate.pkl \
#     --simulation-dt $SIM_DT \
#     --n-days $N_DAYS \
#     --output-basepath $OUTPUT_BASEPATH \
#     --environment-name $ENVIRONMENT_NAME \
#     --code-basepath $MDSINE2_PAPER_CODE_PATH \
#     --queue $QUEUE \
#     --memory $MEM \
#     --n-cpus $N_CPUS \
#     --lsf-basepath $LSF_BASEPATH

# python scripts/run_forward_sim_for_fold.py \
#     --chain output/mdsine2/cv/healthy-cv4/numpy_trace \
#     --validation processed_data/cv/healthy-cv4-validate.pkl \
#     --validation-curr-path ../../processed_data/cv/healthy-cv4-validate.pkl \
#     --simulation-dt $SIM_DT \
#     --n-days $N_DAYS \
#     --output-basepath $OUTPUT_BASEPATH \
#     --environment-name $ENVIRONMENT_NAME \
#     --code-basepath $MDSINE2_PAPER_CODE_PATH \
#     --queue $QUEUE \
#     --memory $MEM \
#     --n-cpus $N_CPUS \
#     --lsf-basepath $LSF_BASEPATH

# python scripts/run_forward_sim_for_fold.py \
#     --chain output/mdsine2/cv/healthy-cv5/numpy_trace \
#     --validation processed_data/cv/healthy-cv5-validate.pkl \
#     --validation-curr-path ../../processed_data/cv/healthy-cv5-validate.pkl \
#     --simulation-dt $SIM_DT \
#     --n-days $N_DAYS \
#     --output-basepath $OUTPUT_BASEPATH \
#     --environment-name $ENVIRONMENT_NAME \
#     --code-basepath $MDSINE2_PAPER_CODE_PATH \
#     --queue $QUEUE \
#     --memory $MEM \
#     --n-cpus $N_CPUS \
#     --lsf-basepath $LSF_BASEPATH

# # UC Cohort
# # ---------
# python scripts/run_forward_sim_for_fold.py \
#     --chain output/mdsine2/cv/uc-cv6/numpy_trace \
#     --validation processed_data/cv/uc-cv6-validate.pkl \
#     --validation-curr-path ../../processed_data/cv/uc-cv6-validate.pkl \
#     --simulation-dt $SIM_DT \
#     --n-days $N_DAYS \
#     --output-basepath $OUTPUT_BASEPATH \
#     --environment-name $ENVIRONMENT_NAME \
#     --code-basepath $MDSINE2_PAPER_CODE_PATH \
#     --queue $QUEUE \
#     --memory $MEM \
#     --n-cpus $N_CPUS \
#     --lsf-basepath $LSF_BASEPATH

# python scripts/run_forward_sim_for_fold.py \
#     --chain output/mdsine2/cv/uc-cv7/numpy_trace \
#     --validation processed_data/cv/uc-cv7-validate.pkl \
#     --validation-curr-path ../../processed_data/cv/uc-cv7-validate.pkl \
#     --simulation-dt $SIM_DT \
#     --n-days $N_DAYS \
#     --output-basepath $OUTPUT_BASEPATH \
#     --environment-name $ENVIRONMENT_NAME \
#     --code-basepath $MDSINE2_PAPER_CODE_PATH \
#     --queue $QUEUE \
#     --memory $MEM \
#     --n-cpus $N_CPUS \
#     --lsf-basepath $LSF_BASEPATH

# python scripts/run_forward_sim_for_fold.py \
#     --chain output/mdsine2/cv/uc-cv8/numpy_trace \
#     --validation processed_data/cv/uc-cv8-validate.pkl \
#     --validation-curr-path ../../processed_data/cv/uc-cv8-validate.pkl \
#     --simulation-dt $SIM_DT \
#     --n-days $N_DAYS \
#     --output-basepath $OUTPUT_BASEPATH \
#     --environment-name $ENVIRONMENT_NAME \
#     --code-basepath $MDSINE2_PAPER_CODE_PATH \
#     --queue $QUEUE \
#     --memory $MEM \
#     --n-cpus $N_CPUS \
#     --lsf-basepath $LSF_BASEPATH

# python scripts/run_forward_sim_for_fold.py \
#     --chain output/mdsine2/cv/uc-cv9/numpy_trace \
#     --validation processed_data/cv/uc-cv9-validate.pkl \
#     --validation-curr-path ../../processed_data/cv/uc-cv9-validate.pkl \
#     --simulation-dt $SIM_DT \
#     --n-days $N_DAYS \
#     --output-basepath $OUTPUT_BASEPATH \
#     --environment-name $ENVIRONMENT_NAME \
#     --code-basepath $MDSINE2_PAPER_CODE_PATH \
#     --queue $QUEUE \
#     --memory $MEM \
#     --n-cpus $N_CPUS \
#     --lsf-basepath $LSF_BASEPATH

# python scripts/run_forward_sim_for_fold.py \
#     --chain output/mdsine2/cv/uc-cv10/numpy_trace \
#     --validation processed_data/cv/uc-cv10-validate.pkl \
#     --validation-curr-path ../../processed_data/cv/uc-cv10-validate.pkl \
#     --simulation-dt $SIM_DT \
#     --n-days $N_DAYS \
#     --output-basepath $OUTPUT_BASEPATH \
#     --environment-name $ENVIRONMENT_NAME \
#     --code-basepath $MDSINE2_PAPER_CODE_PATH \
#     --queue $QUEUE \
#     --memory $MEM \
#     --n-cpus $N_CPUS \
#     --lsf-basepath $LSF_BASEPATH