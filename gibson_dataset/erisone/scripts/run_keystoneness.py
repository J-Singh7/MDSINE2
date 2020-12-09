'''Run keystoneness given the chain.
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
python keystoneness.py \
    --input {chain} \
    --study {study} \
    --leave-out-table {leaveouttable} \
    --leave-out-index {leaveoutindex} \
    --forward-simulate {forward_sim} \
    --make-table {maketable} \
    --sep {sep}
    --simulation-dt {sim_dt} \
    --n-days {n_days} \
    --output-basepath {basepath} \
'''

import mdsine2 as md2
import argparse
import os

if __name__ == '__main__':
    md2.LoggingConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--chain', '-c', type=str, dest='chain',
        help='This is the path of the chain for inference or the folder that contains ' \
             'numpy arrays of the traces for the different parameters')
    parser.add_argument('--study', type=str, dest='study',
        help='Study object to use for initial conditions')
    parser.add_argument('--simulation-dt', type=float, dest='simulation_dt',
        help='Timesteps we go in during forward simulation', default=0.01)
    parser.add_argument('--n-days', type=str, dest='n_days',
        help='Number of days to simulate for', default=None)
    parser.add_argument('--output-basepath', '-o', type=str, dest='basepath',
        help='This is where you are saving the posterior renderings')
    parser.add_argument('--leave-out-table', type=str, dest='leave_out_table',
        help='Table of which taxas to leave out')
    parser.add_argument('--sep', type=str, dest='sep', default=',',
        help='separator for the leave out table')

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

    md2.config.LoggingConfig(level=logging.INFO)
    curr_path_table = os.path.join('../../', args.leave_out_table)
    f = open(curr_path_table, 'r')
    tbl = f.read()
    f.close()
    nlines = len(tbl.split('\n'))

    curr_path_study = os.path.join('../../', args.study)
    study = md2.Study.load(curr_path_table)
    leave_outs = ['none'] + [str(i) for i in range(nlines)]

    lsfdir = args.lsf_basepath
    os.makedirs(lsfdir, exist_ok=True)
    script_path = os.path.join(lsfdir, 'scripts')
    stdout_loc = os.path.join(lsfdir, 'stdout')
    stderr_loc = os.path.join(lsfdir, 'stderr')
    os.makedirs(script_path, exist_ok=True)
    os.makedirs(stdout_loc, exist_ok=True)
    os.makedirs(stderr_loc, exist_ok=True)

    # Dispatch keystoneness
    for leave_out in leave_outs:
        jobname = study.name + '-keystone-{}'.format(leave_out)

        stdout_name = os.path.join(stdout_loc, jobname + '.out')
        stderr_name = os.path.join(stderr_loc, jobname + '.err')
        lsfname = os.path.join(script_path, jobname + '.lsf')

        f = open(lsfname, 'w')
        f.write(lsfstr.format(
            jobname=jobname, stdout_loc=stdout_name, 
            stderr_loc=stderr_name, queue=args.queue, cpus=args.cpus, mem=args.memory,
            environment_name=args.environment_name, code_basepath=args.code_basepath,
            chain=args.chain, study=args.study, leaveouttable=args.leave_out_table,
            leaveoutindex=leave_out, forward_sim=1, maketable=0,
            sep=args.sep, sim_dt=args.simulation_dt, n_days=args.n_days, 
            basepath=args.basepath))
        f.close()
        command = 'bsub < {}'.format(lsfname)
        print(command)
        os.system(command)






