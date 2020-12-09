# ErisOne documentation [Internal use]

### WARNING
__These are for internal use only__. If you try to run these scripts, the job will fail. These jobs run the same as the jobs in the parent folder - they just are submitting as jobs using LSF.

## Programatically generate LSF files for running jobs
Each bash script runs a Python script located in `scripts` that generates a LSF file with the desired parameters.

### Preprocessing
Preprocessing and learning the negative binomial dispersion parameters do not take that much compute power, so we do not need to dispatch these jobs. Instead, see steps in `MDSINE2_Paper/gibson_dataset` for preprocessing (everything before learning the MDSINE2 model). 

### Forward simulation and cross-validation
To run cross validation and forward simulation, run
```bash
./submit_cv_and_tla.sh
```
To only run forward simulation, run
```bash
./submit_tla.sh
```
Computing the errors from a forward simulation does not take that much compute power and can be done on a single interactive node. Use the script `MDSINE2_paper/gibson_dataset/compute_errors_tla.sh`.

### Run the model
Run the model and fixed clustering with two different seeds. Each job first runs the unfixed clustering job then immediately after runs the fixed clustering job. These scripts also visualize the posterior of the jobs.
```bash
./submit_mdsine2.sh
```

### Running on ErisOne/ErisTwo
If running the scripts above you get the following error:
```bash
/bin/bash^M: bad interpreter: No such file or directory
```
this does not mean there is something wrong with the script, this is becuase the permissions are
wrong and/or the script is in the wrong format. This can occur if you work on a windows device and push these scripts to the github. To fix this, run the following commands on the script in question:
```bash
chmod +x script.sh
dos2unix -k -o script.sh
```
`dos2unix` is a standard command on Linux systems so you do not need to install anything.