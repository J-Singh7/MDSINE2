'''Run time lookahead prediction

Author: David Kaplan
Date: 12/01/20
MDSINE2 version: 4.0.6

Input format
------------
There are two different input formats that you can pass in:
    1) MDSINE2.BaseMCMC pickle
    2) Folder of numpy arrays
For most users, passing in the MDSINE2.BaseMCMC file is the way to go. Users
might want to pass in the folder only if they are running jobs in parallel
for this simulation and have many different jobs accessing the data at once. 

You can load the MCMC chain from either a `mcmc.pkl` file or from a folder.
If you load it from a folder, then it must have the following structure
folder/
    growth.npy # np.ndarray(n_gibbs, n_taxa)
    interactions.npy # np.ndarray(n_gibbs, n_taxa, n_taxa)
    perturbations.pkl # Optional, dictionary
        (name of perturbation) str -> dict
            'values' -> np.ndarray (n_gibbs, n_taxa)

Forward simulation
------------------
The default simulation is a full time prediction where we start from the first timepoint
and simulate until the last timepoint. Another thing you can do with this is only 
forward simulate from a subjset of the times. Specify the start time with 
`--start` and the number of days to forward simulate with with `--n-days`. If you
additionally want to save all of the intermediate times within start and end, include
the flag `--save-intermediate-times`.
'''

import argparse
import mdsine2 as md2
from mdsine2.names import STRNAMES
import os
import logging
import numpy as np
import pickle
import time

def forward_simulate(growth, interactions, perturbations, 
    dt, subject, start, n_days, limit_of_detection, full_pred, studyname, 
    basepath, sim_max=None, save_intermediate_times=False):
    '''Forward simulate from day `start` for `n_days` days with data from subject
    `subject`. Record the predicted trajectory (for ever gibb step) in the path
    `predpath` and the ground truth (of the data) in `truthpath`.

    Parameters
    ----------
    growth : np.ndarray(n_gibbs, n_taxa)
        Growth values
    interactons : np.ndarray(n_gibbs, n_taxa, n_taxa)
        Off-diagonal interaction values
    perturbations : dict, None
        None if there are no perturbations
    dt : float
        This is the step size during forward simulation
    subject : md2.Subject
        This is the Subject object that we are getting the data from
    start : float
        This is the day that we start from
    n_days : float
        This is the maximum day to forward simulate to
    limit_of_detection : float
        This is the minimum value the trajectory is set to if the intial 
        conditions are zero
    sim_max : float
        This is the clip of the maximum value during forward simulation
    save_intermediate_times : bool
        If true, We record the times for every timepoint up to number of days for 
        flexibility during postprocessing. Example
            subject.times = [0, 0.5, 1, 1.5, 3, 4, 5]
            start = 0.5
            n_days = 3.5

            We will save the time lookahead for 
                start=0.05, n_days=0.5
                start=0.05, n_days=1.0
                start=0.05, n_days=2.5
                start=0.05, n_days=3.5
            Not just 
                start=0.05, n_days=3.5
    full_pred : bool
        If True, we are doing a full time prediction
    studyname : str
        Name of the study
    basepath : str
        Basepath to save in
    '''
    times = subject.times
    M = subject.matrix()['abs']

    # Get the times and data within the time frame specified
    # ------------------------------------------------------
    if full_pred:
        startidx = 0
        endidx = len(times)
        start = times[0]
        end = times[-1]
        n_days = end - start
    else:
        if start not in times:
            raise ValueError('Start time ({}) not in subject times ({}). This should not ' \
                'happen'.format(start, times))
        
        startidx = np.searchsorted(times, start)
        if startidx == len(times)-1:
            logging.warning('Start time is the last timepoint - nothing to forward simulate.')
            return

        end = start + n_days
        if end not in times:
            if not save_intermediate_times:
                logging.warning('We are not saving intermediate times and the end point ({}) is not ' \
                    'contained in the subject times ({}). Ending'.format(end, times))
                return

            if end > times[-1]:
                end = times[-1]
                endidx = len(times)
            else:

                endidx = None
                for tidx, t in enumerate(times):
                    if t > end:
                        break
                    endidx = tidx+1
                if endidx-1 == startidx:
                    logging.info('`n_days` {} is not large enough from start {} in subject {} ({}). ' \
                        'Ending'.format(n_days, start, subject.name, subject.times))
                    return
                end = times[endidx-1]
            n_days = end - start

        else:
            endidx = np.searchsorted(times, end) + 1 # Add 1 so it is inclusive
        times = times[startidx:endidx]
        M = M[:, startidx:endidx]

    initial_conditions = M[:, 0]
    if np.any(initial_conditions == 0):
        logging.info('{} taxas have a 0 abundance at time {}. Setting to {}'.format(
            np.sum(initial_conditions == 0), start, limit_of_detection))
        initial_conditions[initial_conditions == 0] = limit_of_detection
    initial_conditions = initial_conditions.reshape(-1,1)

    # Make the objects
    # ----------------
    if perturbations is not None:
        perts = []
        pert_starts = []
        pert_ends = []
        for pertname in perturbations:
            perts.append(perturbations[pertname]['value'])
            pert_starts.append(perturbations[pertname]['start'])
            pert_ends.append(perturbations[pertname]['end'])

        # Sort the perturabtions in start order
        idxs = np.argsort(pert_starts)
        tmp_p = []
        tmp_ps = []
        tmp_pe = []
        for idx in idxs:
            tmp_p.append(perts[idx])
            tmp_ps.append(pert_starts[idx])
            tmp_pe.append(pert_ends[idx])
        perts = tmp_p
        pert_starts = tmp_ps
        pert_ends = tmp_pe

    else:
        perts = None
        pert_starts = None
        pert_ends = None
    dyn = md2.model.gLVDynamicsSingleClustering(growth=None, interactions=None, 
        start_day=start, sim_max=sim_max, perturbation_starts=pert_starts,
        perturbation_ends=pert_ends)

    # Forward simulate
    # ----------------
    pred_matrix = np.zeros(shape=(growth.shape[0], growth.shape[1], len(times)))
    start_time = time.time()
    for gibbstep in range(growth.shape[0]):

        if gibbstep % 5 == 0 and gibbstep > 0:
            logging.info('{}/{} - {}'.format(gibbstep,growth.shape[0], 
                time.time()-start_time))
            start_time = time.time()

        dyn.growth = growth[gibbstep]
        dyn.interactions = interactions[gibbstep]
        if perts is not None:
            dyn.perturbations = [pert[gibbstep] for pert in perts]
        
        x = md2.integrate(dynamics=dyn, initial_conditions=initial_conditions,
            dt=dt, n_days=times[-1]+dt, subsample=True, times=times)
        pred_matrix[gibbstep] = x['X']

    # Save the output
    # ---------------
    if full_pred:
        fname = os.path.join(basepath, '{studyname}-{subjname}-full.npy'.format(
            studyname=studyname, subjname=subject.name))
        fname_truth = fname.replace('.npy', '-truth.npy')
        fname_times = fname.replace('.npy', '-times.npy')
        np.save(fname, pred_matrix)
        np.save(fname_truth, M)
        np.save(fname_times, times)
    else:

        if save_intermediate_times:
            # Save all of the sub-times
            for endidx in range(1, len(times)):
                end = times[endidx]
                n_days = end-start
                sliceidx = endidx+1

                truth = M[:, :sliceidx]
                pred = pred_matrix[:, :, :sliceidx]
                ts = times[:sliceidx]

                fname = os.path.join(
                    basepath, '{studyname}-{subjname}-start{start}-ndays{ndays}.npy'.format(
                    studyname=studyname, subjname=subject.name, start=start, ndays=n_days))
                fname_truth = fname.replace('.npy', '-truth.npy')
                fname_times = fname.replace('.npy', '-times.npy')
                np.save(fname, pred)
                np.save(fname_truth, truth)
                np.save(fname_times, ts)

        else:
            fname = os.path.join(
                basepath, '{studyname}-{subjname}-start{start}-ndays{ndays}.npy'.format(
                studyname=studyname, subjname=subject.name, start=start, ndays=n_days))
            fname_truth = fname.replace('.npy', '-truth.npy')
            fname_times = fname.replace('.npy', '-times.npy')
            np.save(fname, pred_matrix)
            np.save(fname_truth, M)
            np.save(fname_times, times)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--input', type=str, dest='input',
        help='Location of input (either folder of the numpy arrays or ' \
            'MDSINE2.BaseMCMC chain)')
    parser.add_argument('--validation', type=str, dest='validation',
        help='Data to do inference with')
    parser.add_argument('--simulation-dt', type=float, dest='simulation_dt',
        help='Timesteps we go in during forward simulation', default=0.01)
    parser.add_argument('--start', type=str, dest='start',
        help='Day to start on', default=None)
    parser.add_argument('--n-days', type=str, dest='n_days',
        help='Number of days to simulate for', default=None)
    parser.add_argument('--limit-of-detection', dest='limit_of_detection',
        help='If any of the taxas have a 0 abundance at the start, then we ' \
            'set it to this value.',default=1e5)
    parser.add_argument('--sim-max', dest='sim_max',
        help='Maximum value', default=1e20)
    parser.add_argument('--output-basepath', '-o', type=str, dest='basepath',
        help='This is where you are saving the posterior renderings')
    parser.add_argument('--save-intermediate-times', type=int, default=0,
        dest='save_intermediate_times', help='If 1, save all of the times between ' \
        'start and start + n_days, not just start + n_days. This is efficient if you ' \
        'are doing many time look ahead predictions at various timepoints')
    
    args = parser.parse_args()
    md2.config.LoggingConfig(level=logging.INFO)
    study = md2.Study.load(args.validation)
    save_intermediate_times = bool(args.save_intermediate_times)
    os.makedirs(args.basepath, exist_ok=True)

    # Set the start and end times
    start_time = args.start
    if start_time is not None:
        if start_time.lower() == 'none':
            start_time = None
        else:
            start_time = float(start_time)

    n_days_total = args.n_days
    if n_days_total is not None:
        if n_days_total.lower() == 'none':
            n_days_total = None
        else:
            n_days_total = float(n_days_total)

    # Get the traces of the parameters
    # --------------------------------
    if '.pkl' in args.input:
        # This is the chain
        logging.info('Input is an MDSINE2.BaseMCMC object')
        mcmc = md2.BaseMCMC.load(args.input)

        growth = mcmc.graph[STRNAMES.GROWTH_VALUE].get_trace_from_disk()
        self_interactions = mcmc.graph[STRNAMES.SELF_INTERACTION_VALUE].get_trace_from_disk()
        interactions = mcmc.graph[STRNAMES.INTERACTIONS_OBJ].get_trace_from_disk()
        interactions[np.isnan(interactions)] = 0
        self_interactions = -np.absolute(self_interactions)
        for i in range(self_interactions.shape[1]):
            interactions[:,i,i] = self_interactions[:, i]

        if mcmc.graph.perturbations is not None:
            logging.info('Perturbations exist')
            perturbations = {}
            for pert in mcmc.graph.perturbations:
                perturbations[pert.name] = {}
                perturbations[pert.name]['value'] = pert.get_trace_from_disk()
                perturbations[pert.name]['value'][np.isnan(perturbations[pert.name]['value'])] = 0

        else:
            logging.info('Did not find perturbations')
            perturbations = None

    else:
        # This is a folder
        logging.info('input is a folder')
        growth = np.load(os.path.join(args.input, 'growth.npy'))
        interactions = np.load(os.path.join(args.input, 'interactions.npy'))
        if os.path.isfile(os.path.join(args.input, 'perturbations.pkl')):
            logging.info('perturbations exist')
            with open(os.path.join(args.input, 'perturbations.pkl'), 'rb') as handle:
                perturbations = pickle.load(handle)

        else:
            logging.info('Did not find perturbations')
            perturbations = None

    if start_time is None and n_days_total is None and not save_intermediate_times:
        full_pred = True
    else:
        full_pred = False

    # Run time lookahead
    for subj in study:
        logging.info('Subject {}'.format(subj.name))
        logging.info('{}'.format(subj.times))
        if start_time is None:
            start = subj.times[0]
        else:
            start = start_time
        if start not in subj.times:
            logging.warning('Start time {} not contained in subject {} times ({}). skipping'.format(
                start, subj.name, subj.times))
            continue
        if n_days_total is None:
            n_days = subj.times[-1] - start
        else:
            n_days = n_days_total

        if perturbations is not None:
            for pert in subj.perturbations:
                perturbations[pert.name]['start'] = pert.starts[subj.name]
                perturbations[pert.name]['end'] = pert.ends[subj.name]

        forward_simulate(
            growth=growth, interactions=interactions,
            perturbations=perturbations, dt=args.simulation_dt, subject=subj,full_pred=full_pred,
            start=start, n_days=n_days, limit_of_detection=args.limit_of_detection, 
            sim_max=args.sim_max, save_intermediate_times=save_intermediate_times,
            studyname=study.name, basepath=args.basepath)
