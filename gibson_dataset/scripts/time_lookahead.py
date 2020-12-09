'''Perform time lookahead for each subject - both full and max-tla
'''

import mdsine2 as md2
import argparse
import os
import numpy as np
import logging

script_format = 'python ../forward_sim.py --input {chain} ' \
    '--validation {validation} ' \
    '--simulation-dt {sim_dt} ' \
    '--start {start} ' \
    '--n-days {n_days} ' \
    '--output-basepath {basepath} ' \
    '--save-intermediate-times {save_inter_times}'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--chain', '-c', type=str, dest='chain',
        help='This is the path of the chain for inference or the folder that contains ' \
             'numpy arrays of the traces for the different parameters')
    parser.add_argument('--validation', type=str, dest='validation',
        help='Data to do inference with')
    parser.add_argument('--simulation-dt', type=float, dest='simulation_dt',
        help='Timesteps we go in during forward simulation', default=0.01)
    parser.add_argument('--n-days', type=str, dest='n_days',
        help='Number of days to simulate for', default=None)
    parser.add_argument('--limit-of-detection', dest='limit_of_detection',
        help='If any of the taxas have a 0 abundance at the start, then we ' \
            'set it to this value.',default=1e5)
    parser.add_argument('--output-basepath', '-o', type=str, dest='basepath',
        help='This is where you are saving the posterior renderings')
    args = parser.parse_args()
    
    # Make the releveant folders and get items
    basepath = args.basepath
    os.makedirs(basepath, exist_ok=True)
    md2.LoggingConfig(level=logging.INFO)

    # Get all the union timepoints within this study object
    study = md2.Study.load(args.validation)
    times = []
    for subj in study:
        times.append(subj.times)
    times = np.sort(np.unique(times))

    # Do a complete lookahead
    command = script_format.format(
        chain=args.chain, validation=args.validation, 
        sim_dt=args.simulation_dt, start=None, n_days=None,
        basepath=basepath, save_inter_times=0)
    logging.info(command)
    os.system(command)

    # Do time lookahead (do not include last time point)
    for start in times[:-1]:
        command = script_format.format(
            chain=args.chain, validation=args.validation, 
            sim_dt=args.simulation_dt, start=start, n_days=args.n_days,
            basepath=basepath, save_inter_times=1)
        logging.info(command)
        os.system(command)

    