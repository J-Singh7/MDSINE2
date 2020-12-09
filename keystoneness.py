'''Compute the keystoneness

Forward simulate the dynamics with taxa(s) left out of the system to see how the steady
state changes

Author: David Kaplan
Date: 12/07/20
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

NOTE: perturbations can be specified in the folder but since they are not used
they are ignored

Leave-out table
---------------
The table passed into `--leave-out-table` is a table that indicates the index of
taxa(s) that should be left out for each keystoneness. Each line represents the taxa(s)
that should be left out. Each line is a separated list of the indexes that 
should be left out. Example:
table.csv
    0,8,4
    1,2,3
    6
    5,9
specify the separator for each element with the `--sep` command. The default separator
is ','

To specify which line you want to compute, set the index with the command `--leave-out-index`.
If you want to compute the steady state without leaving any out, do not specify the 
`--leave-out-index` command. If you want to compute all of the indices, 
use the argument `--leave-out-index all`.

Forward simulate vs. making the table
-------------------------------------
Once you have forward simulated each index of the table, make the table of the steady states.
Taxas that are left out are set to nan.

Default is to forward simulate all of the indexes and to make the table.
To forward simulate, set `--forward-sim 1`. To not forward simulate, set `--forward-sim 0`.
To make the table, set `--make-table 1`. To not make the table, set `--make-table 0`.
'''
import mdsine2 as md2
from mdsine2.names import STRNAMES
import argparse
import logging
import numpy as np
import os
import time
import re
import pandas as pd

def _forward_sim(growth, interactions, initial_conditions, dt, sim_max, n_days):
    '''Forward simulate with the given dynamics
    '''

    dyn = md2.model.gLVDynamicsSingleClustering(growth=None, interactions=None, 
        start_day=0, sim_max=sim_max)
    times = np.arange(n_days+1)
    initial_conditions = initial_conditions.reshape(-1,1)

    start_time = time.time()
    pred_matrix = np.zeros(shape=(growth.shape[0], growth.shape[1], len(times)))
    for gibb in range(growth.shape[0]):
        if gibb % 5 == 0 and gibb > 0:
            logging.info('{}/{} - {}'.format(gibb,growth.shape[0], 
                time.time()-start_time))
            start_time = time.time()
        dyn.growth = growth[gibb]
        dyn.interactions = interactions[gibb]

        x = md2.integrate(dynamics=dyn, initial_conditions=initial_conditions, 
            dt=dt, n_days=times[-1]+dt, subsample=True, times=times)
        pred_matrix[gibb] = x['X']
    return pred_matrix

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--input', type=str, dest='input',
        help='Location of input (either folder of the numpy arrays or ' \
            'MDSINE2.BaseMCMC chain)')
    parser.add_argument('--study', type=str, dest='study',
        help='Study object to use for initial conditions')
    parser.add_argument('--leave-out-table', type=str, dest='leave_out_table',
        help='Table of which taxas to leave out')
    parser.add_argument('--leave-out-index', type=str, dest='leave_out_index',
        help='Index of the table to leave out. If you want to compute all of ' \
             'the lines in `--leave-out-table`, then set to "all". If nothing ' \
             'is passed in, then it will compute the keystoneness with nothing ' \
             'left out.', default='all')
    parser.add_argument('--forward-simulate', '-fs', type=int, dest='forward_simulate',
        help='If 1, compute the forward simulation of the index. If 0, do not', default=1)
    parser.add_argument('--make-table', '-mt', type=int, dest='make_table',
        help='If 1, make the table of the steady states. If 0, do not', default=1)
    parser.add_argument('--sep', type=str, dest='sep', default=',',
        help='separator for the leave out table')
    parser.add_argument('--simulation-dt', type=float, dest='simulation_dt',
        help='Timesteps we go in during forward simulation', default=0.01)
    parser.add_argument('--n-days', type=int, dest='n_days',
        help='Number of days to simulate for', default=60)
    parser.add_argument('--limit-of-detection', dest='limit_of_detection',
        help='If any of the taxas have a 0 abundance at the start, then we ' \
            'set it to this value.',default=1e5, type=float)
    parser.add_argument('--sim-max', dest='sim_max', type=float,
        help='Maximum value', default=1e20)
    parser.add_argument('--output-basepath', '-o', type=str, dest='basepath',
        help='Where to save the output')
    args = parser.parse_args()

    md2.config.LoggingConfig(level=logging.INFO)
    study = md2.Study.load(args.study)
    basepath = os.path.join(args.basepath, study.name)
    os.makedirs(basepath, exist_ok=True)

    # Get the initial conditions
    # --------------------------
    M = study.matrix(dtype='abs', agg='mean', times='intersection')
    initial_conditions_master = M[:,0]
    initial_conditions_master[initial_conditions_master==0] = args.limit_of_detection

    # Load table
    # ----------
    f = open(args.leave_out_table, 'r')
    tbl = f.read()
    f.close()
    lines = tbl.split('\n')

    # Remove empty lines
    _lines = []
    for line in lines:
        if line == '':
            continue
        _lines.append(line)
    lines = _lines

    # Index out if necessary
    # ----------------------
    leave_out_index = args.leave_out_index
    if md2.isstr(leave_out_index):
        leave_out_index = leave_out_index.lower()
    if leave_out_index is None or leave_out_index == 'none':
        logging.info('No index provided, not indexing out any taxa')
        idxs = [None]
    elif leave_out_index == 'all':
        logging.info('"all" provided as index. Iterating over each index')
        idxs = np.arange(len(lines))
        idxs = [None] + idxs.tolist()
    else:
        try:
            idxs = int(leave_out_index)
        except:
            logging.critical('--leave-out-index ({}) not recognized as an index'.format(
                leave_out_index))
            raise
        idxs = [idxs]

    for idx in idxs:
        if idx is not None:
            if idx >= len(lines):
                raise ValueError('index {} is out range ({} total)'.format(idx, len(lines)))

    if args.forward_simulate == 1:
        logging.info('Forward simulating')
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
            growth_master = growth
            interactions_master = interactions

        else:
            # This is a folder
            logging.info('input is a folder')
            growth_master = np.load(os.path.join(args.input, 'growth.npy'))
            interactions_master = np.load(os.path.join(args.input, 'interactions.npy'))

        # Forward simulate
        for idx in idxs:
            if idx is not None:
                oidxs = [int(ele) for ele in lines[idx].split(args.sep)]
                logging.info('indexing out {}'.format(oidxs))

                mask = np.ones(len(initial_conditions_master), dtype=bool)
                mask[oidxs] = False
                growth = growth_master[:, mask]
                interactions = interactions_master[:, mask, :]
                interactions = interactions[:, :, mask]
                initial_conditions = initial_conditions_master[mask]
            else:
                logging.info('not indexing out anything')
                growth = growth_master
                interactions = interactions_master
                initial_conditions = initial_conditions_master

            pred_matrix = _forward_sim(growth=growth, interactions=interactions, 
                initial_conditions=initial_conditions, dt=args.simulation_dt,
                sim_max=args.sim_max, n_days=args.n_days)

            # Save the forward sims
            if idx is None:
                name = 'none'
            else:
                name = str(idx)
            fname = os.path.join(basepath, 'study{}-lo{}-forward-sims.npy'.format(study.name, name))
            np.save(fname, pred_matrix)
    if args.make_table == 1:
        # Make the table
        # --------------
        logging.info('Make the table')

        re_find = re.compile(r'^study(.*)-lo(.*)-forward-sims.npy$')

        fnames = os.listdir(basepath)
        data = []
        idxs = []
        for fname in fnames:
            try:
                studyname, leaveout = re_find.findall(fname)[0]
            except:
                logging.warning('{} not recognized as a pattern. skipping'.format(fname))
                continue
            if studyname != study.name:
                logging.warning('{} not the same study as input {}. skipping'.format(
                    studyname, study.name))
                continue

            if leaveout == 'none':
                leaveout = -1
            else:
                leaveout = int(leaveout)

            M = np.load(os.path.join(basepath, fname))
            steady_state = np.nanmean(M[:, :, 0], axis=0)
            
            if leaveout != -1:
                oidxs = [int(ele) for ele in lines[leaveout].split(args.sep)]
                mask = np.ones(len(study.taxas), dtype=bool)
                mask[oidxs] = False
                ss = np.zeros(len(study.taxas)) * np.nan
                ss[mask] = steady_state
                steady_state = ss
            
            data.append(steady_state)
            idxs.append(leaveout)
        
        if len(data) != len(lines) + 1:
            raise ValueError('You are making the table. {} total lines were found instead of {}'.format( 
                len(data), len(lines)+1))
        
        # order the table
        data = np.asarray(data)
        idxs = np.argsort(idxs)
        data = data[idxs, :]

        columns = [taxa.name for taxa in study.taxas]
        index = ['base'] + lines

        df = pd.DataFrame(data, index=index, columns=columns)
        df.to_csv(os.path.join(basepath, 'table.tsv'), sep='\t', index=True, header=True)