'''Make the lsf file for cross validation.

WARNING: THIS FILE ONLY WORKS IF THE OS CONTAINS AN LSF JOB SUBMISSION SYSTEM.
THIS IS A INTERNAL DOCUMENT. FOR IDENTICAL RESULTS THAT DO NOT REQUIRE RUNNING 
LSF, RUN THE SCRIPT `MDSINE2/gibson_dataset/run_cv.sh`.
'''

lsfstr = '''#!/bin/bash
#BSUB -J {jobname}
#BSUB -o {stdout_loc}
#BSUB -e {stderr_loc}

#BSUB -q {cv_queue}
#BSUB -n {cv_cpus}
#BSUB -M {cv_mem}
#BSUB -R rusage[mem={cv_mem}]

# 
echo '---PROCESS RESOURCE LIMITS---'
ulimit -a
echo '---SHARED LIBRARY PATH---'
echo $LD_LIBRARY_PATH
echo '---APPLICATION SEARCH PATH:---'
echo $PATH
echo '---LSF Parameters:---'
printenv | grep '^LSF'
echo '---LSB Parameters:---'
printenv | grep '^LSB'
echo '---LOADED MODULES:---'
module list
echo '---SHELL:---'
echo $SHELL
echo '---HOSTNAME:---'
hostname
echo '---GROUP MEMBERSHIP (files are created in the first group listed):---'
groups
echo '---DEFAULT FILE PERMISSIONS (UMASK):---'
umask
echo '---CURRENT WORKING DIRECTORY:---'
pwd
echo '---DISK SPACE QUOTA---'
df .
echo '---TEMPORARY SCRATCH FOLDER ($TMPDIR):---'
echo $TMPDIR


# Load the environment
module load anaconda/default
source activate {environment_name}
cd {code_basepath}

# Run a fold in cross validation
python run_cross_validation.py \
    --dataset {dset_fileloc} \
    --cv-basepath {cv_basepath} \
    --dset-basepath {dset_basepath} \
    --negbin {negbin_run} \
    --seed {seed} \
    --burnin {burnin} \
    --n-samples {n_samples} \
    --checkpoint {checkpoint} \
    --multiprocessing {mp} \
    --leave-out-subject {leave_out_subject} \
    --interaction-ind-prior {interaction_prior} \
    --perturbation-ind-prior {perturbation_prior}

# Make the posterior as numpy arrays
python gibson_dataset/scripts/convert_trace_to_numpy.py \
    --chain {chain_path} \
    --output-basepath {numpy_basepath} \
    --section posterior

# Visualize the posterior
python ../../step_6_visualize_mdsine2.py \
    --chain {chain_path} \
    --section posterior \
    --output-basepath {posterior_basepath} \
    --fixed-clustering 0

# Compute forward simulations for this fold
python gibson_dataset/erisone/scripts/run_forward_sim_for_fold.py \
    --chain {numpy_basepath} \
    --validation {validation_subject} \
    --n-days {max_tla} \
    --output-basepath {tla_basepath} \
    --environment-name {environment_name} \
    --code-basepath {code_basepath} \
    --queue {tla_queue} \
    --memory {tla_mem} \
    --n-cpus {tla_cpus} \
    --lsf-basepath {tla_lsf_basepath}
'''

import mdsine2 as md2
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--dataset', '-d', type=str, dest='dataset',
        help='This is the Gibson dataset we want to do cross validation on. This path' \
             ' is relative to the code basepath.')
    parser.add_argument('--dataset-curr-path', '-dcp', type=str, dest='dataset_curr_path',
        help='This is the Gibson dataset we want to do cross validation on. This path ' \
             'is relative to the erisone path')
    parser.add_argument('--cv-basepath', '-o', type=str, dest='output_basepath',
        help='This is the basepath to save the output')
    parser.add_argument('--dset-basepath', '-db', type=str, dest='input_basepath',
        help='This is the basepath to load and save the cv datasets')
    parser.add_argument('--lsf-basepath', '-l', type=str, dest='lsf_basepath',
        help='This is the basepath to save the lsf files', default='lsf_files/')
    parser.add_argument('--leave-out-subject', '-lo', type=str, dest='leave_out_subj',
        help='This is the subject to leave out')
    parser.add_argument('--negbin', type=str, dest='negbin',
        help='This is the MCMC object that was run to learn a0 and a1')
    parser.add_argument('--seed', '-s', type=int, dest='seed',
        help='This is the seed to initialize the inference with')
    parser.add_argument('--burnin', '-nb', type=int, dest='burnin',
        help='How many burn-in Gibb steps for Markov Chain Monte Carlo (MCMC)')
    parser.add_argument('--n-samples', '-ns', type=int, dest='n_samples',
        help='Total number Gibb steps to perform during MCMC inference')
    parser.add_argument('--checkpoint', '-c', type=int, dest='checkpoint',
        help='How often to write the posterior to disk. Note that `--burnin` and ' \
             '`--n-samples` must be a multiple of `--checkpoint` (e.g. checkpoint = 100, ' \
             'n_samples = 600, burnin = 300)')
    parser.add_argument('--multiprocessing', '-mp', type=int, dest='mp',
        help='If 1, run the inference with multiprocessing. Else run on a single process',
        default=0)
    parser.add_argument('--interaction-ind-prior', '-ip', type=str, dest='interaction_prior',
        help='Prior of the indicator of the interactions')
    parser.add_argument('--perturbation-ind-prior', '-pp', type=str, dest='perturbation_prior',
        help='Prior of the indicator of the perturbations')

    # TLA arguments
    parser.add_argument('--max-tla', type=int, dest='max_tla',
        help='Maximum time for time lookahead')
    parser.add_argument('--simulation-dt', type=float, dest='simulation_dt',
        help='Timesteps we go in during forward simulation', default=0.01)
    
    #ErisOne arguments
    parser.add_argument('--environment-name', dest='environment_name', type=str,
        help='Name of the conda environment to activate when the job starts')
    parser.add_argument('--code-basepath', type=str, dest='code_basepath',
        help='Where the `run_cross_validation` script is located')
    parser.add_argument('--cv-queue', type=str, dest='cv_queue',
        help='ErisOne queue this job gets submitted to for cross-validation')
    parser.add_argument('--cv-memory', type=str, dest='cv_memory',
        help='Amount of memory to reserve on ErisOne for cross-validation')
    parser.add_argument('--cv-n-cpus', type=str, dest='cv_cpus',
        help='Number of cpus to reserve on ErisOne for cross-validation')
    parser.add_argument('--tla-queue', type=str, dest='tla_queue',
        help='ErisOne queue this job gets submitted to')
    parser.add_argument('--tla-memory', type=str, dest='tla_memory',
        help='Amount of memory to reserve on ErisOne')
    parser.add_argument('--tla-n-cpus', type=str, dest='tla_cpus',
        help='Number of cpus to reserve on ErisOne')

    args = parser.parse_args()
    study = md2.Study.load(args.dataset_curr_path)
    dset = study.name

    # Make directory for lsf files
    os.makedirs(args.lsf_basepath, exist_ok=True)
    lsf_basepath = os.path.join(args.lsf_basepath, dset)
    os.makedirs(lsf_basepath, exist_ok=True)

    jobname = dset + '-cv' + args.leave_out_subj
    script_path = os.path.join(lsf_basepath, 'scripts')
    stdout_loc = os.path.join(lsf_basepath, 'stdout', jobname + '.out')
    stderr_loc = os.path.join(lsf_basepath, 'stderr', jobname + '.err')
    os.makedirs(script_path, exist_ok=True)
    os.makedirs(stdout_loc, exist_ok=True)
    os.makedirs(stderr_loc, exist_ok=True)

    # Make parameters for time-lookahead
    chain_path = os.path.join(args.output_basepath, jobname, 'mcmc.pkl')
    numpy_basepath = os.path.join(args.output_basepath, jobname, 'numpy_trace')
    posterior_basepath = os.path.join(args.output_basepath, jobname, 'posterior')
    validation_subject = os.path.join(args.input_basepath, jobname + '-validate.pkl')
    tla_basepath = os.path.join(args.output_basepath, 'forward_sims')
    tla_lsf_basepath = os.path.join(lsf_basepath, 'tla')

    os.makedirs(numpy_basepath, exist_ok=True)
    os.makedirs(tla_basepath, exist_ok=True)
    os.makedirs(tla_lsf_basepath, exist_ok=True)

    lsfname = os.path.join(script_path, jobname + '.lsf')
    f = open(lsfname, 'w')
    f.write(lsfstr.format(
        jobname=jobname, stdout_loc=stdout_loc, stderr_loc=stderr_loc,
        environment_name=args.environment_name, 
        cv_queue=args.cv_queue, cv_cpus=args.cv_cpus, cv_mem=args.cv_memory,
        code_basepath=args.code_basepath,
        dset_fileloc=args.dataset, cv_basepath=args.output_basepath,
        dset_basepath=args.input_basepath, negbin_run=args.negbin, 
        seed=args.seed, burnin=args.burnin, n_samples=args.n_samples,
        checkpoint=args.checkpoint, mp=args.mp, 
        leave_out_subject=args.leave_out_subj,
        interaction_prior=args.interaction_prior,
        perturbation_prior=args.perturbation_prior,
        chain_path=chain_path, numpy_basepath=numpy_basepath,
        posterior_basepath=posterior_basepath, validation_subject=validation_subject, 
        max_tla=args.max_tla, tla_basepath=tla_basepath, 
        tla_queue=args.tla_queue, tla_mem=args.tla_memory, 
        tla_cpus=args.tla_cpus, tla_lsf_basepath=tla_lsf_basepath))
    f.close()
    command = 'bsub < {}'.format(lsfname)
    print(command)
    os.system(command)