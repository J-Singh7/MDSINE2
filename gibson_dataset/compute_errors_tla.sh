
OUTPUT_BASEPATH="../output/mdsine2/cv/forward_sims"
ERROR_TABLE_PATH="../output/mdsine2/cv/forward_sims/errors.tsv"

# Compute errors
python ../compute_forward_sim_error.py \
    --input ${OUTPUT_BASEPATH} \
    --output ${ERROR_TABLE_PATH} \
    --error RMSE relRMSE spearman