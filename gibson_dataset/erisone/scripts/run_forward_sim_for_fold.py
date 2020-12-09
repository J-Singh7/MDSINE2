'''Run time lookahead prediction and a full trajectory prediction given the chain 
and data. Submit a job for every forward prediction task.

This is called after cross validation is done for the fold and is automatically
called from the lsf script in `run_cv.py`
'''

lsfstr = '''#!/bin/bash
#BSUB -J {jobname}
#BSUB -o {stdout_loc}
#BSUB -e {stderr_loc}

#BSUB -q {queue}
#BSUB -n {cpus}
#BSUB -M {mem}
#BSUB -R rusage[mem={mem}]

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

# Run time lookahead
python forward_sim.py \
    --input {chain} \
    --validation {validation_path} \
    --simulation-dt {sim_dt} \
    --start {start} \
    --n-days {n_days} \
    --output-basepath {basepath} \
    --save-intermediate-times {save_intermed_times}
'''

import mdsine2 as md2
import argparse
import os
import numpy as np
import logging

if __name__ == '__main__':
    md2.LoggingConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--chain', '-c', type=str, dest='chain',
        help='This is the path of the chain for inference or the folder that contains ' \
             'numpy arrays of the traces for the different parameters')
    parser.add_argument('--validation', type=str, dest='validation',
        help='Data to do inference with. This we use in the lsf file')
    parser.add_argument('--validation-curr-path', type=str, dest='validation_curr_path',
        help='Data with a path relative to erisone folder. This is used in this script',
        default=None)
    parser.add_argument('--simulation-dt', type=float, dest='simulation_dt',
        help='Timesteps we go in during forward simulation', default=0.01)
    parser.add_argument('--n-days', type=str, dest='n_days',
        help='Number of days to simulate for', default=None)
    parser.add_argument('--limit-of-detection', dest='limit_of_detection',
        help='If any of the taxas have a 0 abundance at the start, then we ' \
            'set it to this value.',default=1e5)
    parser.add_argument('--output-basepath', '-o', type=str, dest='basepath',
        help='This is where you are saving the posterior renderings')

    # ErisOne parameters
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
    parser.add_argument('--lsf-basepath', '-l', type=str, dest='lsf_basepath',
        help='This is the basepath to save the lsf files', default='lsf_files/')

    args = parser.parse_args()
    n_days =args.n_days
    basepath = args.basepath

    if args.validation_curr_path is None:
        logging.info('`validation_curr_path` not passed in. St to the same as validation')
        args.validation_curr_path = args.validation

    # Get all the union timepoints within this study object
    study = md2.Study.load(args.validation_curr_path )
    times = []
    for subj in study:
        times.append(subj.times)
    times = np.sort(np.unique(times))

    lsfdir = args.lsf_basepath
    os.makedirs(lsfdir, exist_ok=True)
    script_path = os.path.join(lsfdir, 'scripts')
    stdout_loc = os.path.join(lsfdir, 'stdout')
    stderr_loc = os.path.join(lsfdir, 'stderr')
    os.makedirs(script_path, exist_ok=True)
    os.makedirs(stdout_loc, exist_ok=True)
    os.makedirs(stderr_loc, exist_ok=True)

    # Do time lookahead (do not include last time point)
    for start in times[:-1]:
        jobname = study.name + '-{}-{}'.format(start, n_days)

        stdout_name = os.path.join(stdout_loc, jobname + '.out')
        stderr_name = os.path.join(stderr_loc, jobname + '.err')
        lsfname = os.path.join(script_path, jobname + '.lsf')
        f = open(lsfname, 'w')
        f.write(lsfstr.format(
            jobname=jobname, stdout_loc=stdout_name, stderr_loc=stderr_name,
            queue=args.queue, cpus=args.cpus, mem=args.memory,
            environment_name=args.environment_name, code_basepath=args.code_basepath,
            chain=args.chain, validation_path=args.validation, sim_dt=args.simulation_dt,
            start=start, n_days=n_days, basepath=basepath, save_intermed_times=1))
        f.close()
        command = 'bsub < {}'.format(lsfname)
        print(command)
        os.system(command)

    # Do full simulation
    jobname = study.name + '-full'
    stdout_name = os.path.join(stdout_loc, jobname + '.out')
    stderr_name = os.path.join(stderr_loc, jobname + '.err')
    lsfname = os.path.join(script_path, jobname + '.lsf')
    f = open(lsfname, 'w')
    f.write(lsfstr.format(
        jobname=jobname, stdout_loc=stdout_name, stderr_loc=stderr_name,
        queue=args.queue, cpus=args.cpus, mem=args.memory,
        environment_name=args.environment_name, code_basepath=args.code_basepath,
        chain=args.chain, validation_path=args.validation, sim_dt=args.simulation_dt,
        start=None, n_days=None, basepath=basepath, save_intermed_times=0))
    f.close()
    command = 'bsub < {}'.format(lsfname)
    print(command)
    os.system(command)



