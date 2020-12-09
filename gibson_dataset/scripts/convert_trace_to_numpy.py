'''Convert the parameters in a MCMC chain object into a folder of numpy arrays.
This is used so that you can forward simulate from a trace without having to worry
about permission in HDF5.
'''

import mdsine2 as md2
from mdsine2.names import STRNAMES
import numpy as np
import argparse
import os
import logging
import pickle

if __name__ == '__main__':

    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--chain', '-c', type=str, dest='chain',
        help='This is the chain to save')
    parser.add_argument('--section', '-s', type=str, dest='section',
        help='This is the section of the trace that we are saving',
        default='posterior')
    parser.add_argument('--output-basepath', '--basepath', '-o', type=str, 
        dest='basepath',
        help='This is the folder that we should save the numpy arrays')
    args = parser.parse_args()

    md2.LoggingConfig(level=logging.INFO)
    basepath = args.basepath
    section = args.section
    os.makedirs(basepath, exist_ok=True)

    # Save the different parameters of the chain
    mcmc = md2.BaseMCMC.load(args.chain)

    # growth
    logging.info('Loading growth')
    growth_trace = mcmc.graph[STRNAMES.GROWTH_VALUE].get_trace_from_disk(section=section)
    np.save(os.path.join(basepath, 'growth.npy'), growth_trace)
    growth_trace = None

    # interactions
    logging.info('loading interactions')
    si_trace = -np.absolute(mcmc.graph[STRNAMES.SELF_INTERACTION_VALUE].get_trace_from_disk(section=section))
    interactions_trace = mcmc.graph[STRNAMES.INTERACTIONS_OBJ].get_trace_from_disk(section=section)

    interactions_trace[np.isnan(interactions_trace)] = 0
    for i in range(len(mcmc.graph.data.taxas)):
        interactions_trace[:,i,i] = si_trace[:,i]

    np.save(os.path.join(basepath, 'interactions.npy'), interactions_trace)
    interactions_trace = None
    si_trace = None

    # perturbations
    if mcmc.graph.perturbations is not None:
        logging.info('Loading perturbations')

        perts = {}
        for pert in mcmc.graph.perturbations:
            logging.info('Perturbation {}'.format(pert.name))
            perts[pert.name] = {}
            perts[pert.name]['value'] = pert.get_trace_from_disk()
            perts[pert.name]['value'][np.isnan(perts[pert.name]['value'])] = 0

        pert_fname = os.path.join(basepath, 'perturbations.pkl')
        with open(pert_fname, 'wb') as handle:
            pickle.dump(perts, handle, protocol=pickle.HIGHEST_PROTOCOL)

    else:
        logging.info('There are no perturbations')