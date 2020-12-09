'''Make the lsf files for learning the parameters for the MDSINE2 model.
This first runs the non-fixed clustering and then fixed clustering immediately after.

WARNING: THIS FILE ONLY WORKS IF THE OS CONTAINS AN LSF JOB SUBMISSION SYSTEM.
THIS IS A INTERNAL DOCUMENT. FOR IDENTICAL RESULTS THAT DO NOT REQUIRE RUNNING 
LSF, RUN THE SCRIPT `MDSINE2/gibson_dataset/run_msine2.sh` and 
`MDSINE2/gibson_dataset/run_mdsien2_fixed_clustering`.
'''

lsfstr = '''#!/bin/bash
#BSUB -J {jobname}
#BSUB -o {stdout_loc}
#BSUB -e {stderr_loc}

#BSUB -q {queue}
#BSUB -n {cpus}
#BSUB -M {mem}
#BSUB -R rusage[mem={mem}]

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

# Run the model
# -------------
python step_5_infer_mdsine2.py \
    --input {dset_fileloc} \
    --negbin {negbin_run} \
    --seed {seed} \
    --burnin {burnin} \
    --n-samples {n_samples} \
    --checkpoint {checkpoint} \
    --multiprocessing {mp} \
    --rename-study {rename_study} \
    --output-basepath {basepath} \
    --interaction-ind-prior {interaction_prior} \
    --perturbation-ind-prior {perturbation_prior}
# Plot the posterior
python step_6_visualize_mdsine2.py \
    --chain {chain_loc} \
    --output-basepath {posterior_path} \
    --section posterior \
    --fixed-clustering 0

# Run the model that was just learned with fixed clustering
# ---------------------------------------------------------
python step_5_infer_mdsine2.py \
    --input {dset_fileloc} \
    --fixed-clustering {chain_loc} \
    --negbin {negbin_run} \
    --seed {seed} \
    --burnin {burnin} \
    --n-samples {n_samples} \
    --checkpoint {checkpoint} \
    --multiprocessing {mp} \
    --rename-study {rename_study} \
    --output-basepath {fixed_basepath} \
    --interaction-ind-prior {interaction_prior} \
    --perturbation-ind-prior {perturbation_prior}
# Plot the posterior of fixed clustering
python step_6_visualize_mdsine2.py \
    --chain {fixed_chain_loc} \
    --output-basepath {fixed_posterior_path} \
    --section posterior \
    --fixed-clustering 1
'''

import mdsine2 as md2
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--dataset', '-d', type=str, dest='dataset',
        help='This is the Gibson dataset we are performing inference on')
    parser.add_argument('--lsf-basepath', '-l', type=str, dest='lsf_basepath',
        help='This is the basepath to save the lsf files', default='lsf_files/')
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
    parser.add_argument('--rename-study', type=str, dest='rename_study',
        help='Rename the name of the study to this', default=None)
    parser.add_argument('--output-basepath', type=str, dest='basepath',
        help='Output of the model', default=None)
    parser.add_argument('--fixed-output-basepath', type=str, dest='fixed_basepath',
        help='Output of the fixed-clustering model', default=None)
    parser.add_argument('--interaction-ind-prior', '-ip', type=str, dest='interaction_prior',
        help='Prior of the indicator of the interactions')
    parser.add_argument('--perturbation-ind-prior', '-pp', type=str, dest='perturbation_prior',
        help='Prior of the indicator of the perturbations')
    
    # Erisone Parameters
    parser.add_argument('--environment-name', dest='environment_name', type=str,
        help='Name of the conda environment to activate when the job starts')
    parser.add_argument('--code-basepath', type=str, dest='code_basepath',
        help='Where the `run_cross_validation` script is located')
    parser.add_argument('--queue', '-q', type=str, dest='queue',
        help='ErisOne queue this job gets submitted to')
    parser.add_argument('--memory', '-mem', type=str, dest='memory',
        help='Amount of memory to reserve on ErisOne')
    parser.add_argument('--n-cpus', '-cpus', type=str, dest='cpus',
        help='Number of cpus to reserve on ErisOne')
    
    args = parser.parse_args()

    # Make the arguments
    jobname = args.rename_study
    chain_loc = os.path.join(args.basepath, jobname, 'mcmc.pkl')
    posterior_path = os.path.join(args.basepath, jobname, 'posterior')

    fixed_chain_loc = os.path.join(args.fixed_basepath, jobname, 'mcmc.pkl')
    fixed_posterior_path = os.path.join(args.fixed_basepath, jobname, 'posterior')

    lsfdir = args.lsf_basepath

    script_path = os.path.join(lsfdir, 'scripts')
    stdout_loc = os.path.abspath(os.path.join(lsfdir, 'stdout'))
    stderr_loc = os.path.abspath(os.path.join(lsfdir, 'stderr'))
    os.makedirs(script_path, exist_ok=True)
    os.makedirs(stdout_loc, exist_ok=True)
    os.makedirs(stderr_loc, exist_ok=True)
    stdout_loc = os.path.join(stdout_loc, jobname + '.out')
    stderr_loc = os.path.join(stderr_loc, jobname + '.err')

    os.makedirs(lsfdir, exist_ok=True)
    lsfname = os.path.join(script_path, jobname + '.lsf')

    f = open(lsfname, 'w')
    f.write(lsfstr.format(
        jobname=jobname, stdout_loc=stdout_loc, stderr_loc=stderr_loc, 
        environment_name=args.environment_name,
        code_basepath=args.code_basepath, queue=args.queue, cpus=args.cpus, 
        mem=args.memory, dset_fileloc=args.dataset, 
        negbin_run=args.negbin, seed=args.seed, burnin=args.burnin, 
        n_samples=args.n_samples, checkpoint=args.checkpoint, mp=args.mp,
        rename_study=args.rename_study, basepath=args.basepath, 
        chain_loc=chain_loc, posterior_path=posterior_path, 
        fixed_basepath=args.fixed_basepath, fixed_chain_loc=fixed_chain_loc, 
        fixed_posterior_path=fixed_posterior_path,
        interaction_prior=args.interaction_prior,
        perturbation_prior=args.perturbation_prior))
    f.close()
    command = 'bsub < {}'.format(lsfname)
    print(command)
    os.system(command)
