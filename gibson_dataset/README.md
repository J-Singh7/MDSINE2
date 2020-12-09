# Replication of MDSINE2 results

### Preprocess the data (Note that this has already been done for this dataset)

Preprocessing and agglomeration. Note that when running data from scratch, manual steps are involved (phylogenetic placement, etc.). All of the outputs from preprocessing are already provided here for you. For information on manual steps see [internal_doc_for_manual_steps.md](internal_doc_for_manual_steps.md) before running the following command:
```bash
./preprocessing_agglomeration.sh
```

### Filtering and visualizing the data (the tutorials start here)
Visualize the OTUs and filter the data 
```bash
./plot_aggregates.sh
./plot_phylogenetic_subtrees.sh
./preprocessing_filtering.sh
```

### Learn Negative binomial dispersion parameters
Learn the negative binomial dispersion parameters
```bash
./learn_negbin.sh
```

### Cross-validation and forward simulation
Order of scripts from start to finish of running forward simulation and cross validation:
```bash
./run_cv.sh
./run_tla.sh
./compute_errors_tla.sh
```

### Learning parameters of MDSINE2
Order of scripts from start to finish of generating the posteriors

```bash
./run_mdsine2.sh
./run_mdsine2_fixed_clustering.sh
```

### Post-processing
Once `run_mdsine2.sh` has finished running, you can perform keystoneness and the perturbation analysis
```bash
./compute_keystoneness.sh
./compute_perturbation_analysis.sh
```

### Making figures
Once cross-validation and learning the parameters are done, you can generate the figures used in the paper:
```bash
./make_figures.sh
```