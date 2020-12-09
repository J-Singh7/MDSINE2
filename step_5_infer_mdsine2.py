'''Run MDSINE2 inference

Author: David Kaplan
Date: 11/30/20
MDSINE2 version: 4.0.6

Fixed Clustering
----------------
To run inference with fixed clustering, use the parameter `--fixed-clustering` where
this is the location of the MCMC object ob the inference. This will automatically set the
parameters for the clustering intialization and set learning turned off for the 
cluster assignments.

Priors on the sparsity
----------------------
Set the sparsity of the prior indicator of the interactions and perturbations with the
arguments `--interaction-ind-prior` and `perturbation-ind-prior`.
'''
import argparse
import mdsine2 as md2
import logging
import os
import time
from mdsine2.names import STRNAMES

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--input', '-i', type=str, dest='input',
        help='This is the dataset to do inference with.')
    parser.add_argument('--fixed-clustering', type=str, dest='fixed_clustering',
        help='If you are running fixed clustering, this is the location of the chain ' \
             'that you are setting from.', default=None)
    parser.add_argument('--negbin', type=str, dest='negbin', nargs='+',
        help='If there is a single argument, then this is the MCMC object that was run to ' \
             'learn a0 and a1. If there are two arguments passed, these are the a0 and a1 ' \
             'of the negative binomial dispersion parameters. Example: ' \
             '--negbin /path/to/negbin/mcmc.pkl. Example: ' \
             '--negbin 0.0025 0.025')
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
    parser.add_argument('--basepath', '--output-basepath','-b', type=str, dest='basepath',
        help='This is folder to save the output of inference')
    parser.add_argument('--multiprocessing', '-mp', type=int, dest='mp',
        help='If 1, run the inference with multiprocessing. Else run on a single process',
        default=0)
    parser.add_argument('--rename-study', type=str, dest='rename_study',
        help='Specify the name of the study to set', default=None)
    parser.add_argument('--interaction-ind-prior', '-ip', type=str, dest='interaction_prior',
        help='Prior of the indicator of the interactions')
    parser.add_argument('--perturbation-ind-prior', '-pp', type=str, dest='perturbation_prior',
        help='Prior of the indicator of the perturbations')
    args = parser.parse_args()

    # 1) load dataset
    logging.info('Loading dataset {}'.format(args.input))
    study = md2.Study.load(args.input)
    if args.rename_study is not None:
        if args.rename_study.lower() != 'none':
            study.name = args.rename_study

    # 2) Load the model parameters
    os.makedirs(args.basepath, exist_ok=True)
    basepath = os.path.join(args.basepath, study.name)
    os.makedirs(basepath, exist_ok=True)
    md2.config.LoggingConfig(level=logging.INFO, basepath=basepath)

    # Load the negative binomial parameters
    if len(args.negbin) == 1:
        negbin = md2.BaseMCMC.load(args.negbin[0])
        a0 = md2.summary(negbin.graph[STRNAMES.NEGBIN_A0])['mean']
        a1 = md2.summary(negbin.graph[STRNAMES.NEGBIN_A1])['mean']
        
    elif len(args.negbin) == 2:
        a0 = float(args.negbin[0])
        a1 = float(args.negbin[1])
    else:
        raise ValueError('There must be only 1 or two arguments. use `python ' \
            'step_5_infer_mdsine2.py --help` for options')

    logging.info('Setting a0 = {:.4E}, a1 = {:.4E}'.format(a0,a1))

    # 3) Begin inference
    params = md2.config.MDSINE2ModelConfig(
        basepath=basepath, data_seed=args.seed, init_seed=args.seed, 
        burnin=args.burnin, n_samples=args.n_samples, negbin_a1=a1, 
        negbin_a0=a0, checkpoint=args.checkpoint)
    # Run with multiprocessing if necessary
    if args.mp == 1:
        params.MP_FILTERING = 'full'
        params.MP_CLUSTERING = 'full-4'

    # Change parameters if there is fixed clustering
    if args.fixed_clustering:
        params.LEARN[STRNAMES.CLUSTERING] = False
        params.LEARN[STRNAMES.CONCENTRATION] = False

        params.INITIALIZATION_KWARGS[STRNAMES.CLUSTERING]['value_option'] = 'fixed-clustering'
        params.INITIALIZATION_KWARGS[STRNAMES.CLUSTERING]['value'] = args.fixed_clustering

    # Set the sparsities
    if args.interaction_prior is None:
        raise ValueError('Must specify `--interaction-ind-prior`')
    params.INITIALIZATION_KWARGS[STRNAMES.CLUSTER_INTERACTION_INDICATOR_PROB]['hyperparam_option'] = \
        args.interaction_prior
    if args.perturbation_prior is None:
        raise ValueError('Must specify `--perturbation-ind-prior`')
    params.INITIALIZATION_KWARGS[STRNAMES.PERT_INDICATOR_PROB]['hyperparam_option'] = \
        args.perturbation_prior


    mcmc = md2.initialize_graph(params=params, graph_name=study.name, subjset=study)
    mdata_fname = os.path.join(params.MODEL_PATH, 'metadata.txt')
    params.make_metadata_file(fname=mdata_fname)

    start_time = time.time()
    mcmc = md2.run_graph(mcmc, crash_if_error=True)

    # Record how much time inference took
    t = time.time() - start_time
    t = t / 3600 # Convert to hours

    f = open(mdata_fname, 'a')
    f.write('\n\nTime for inference: {} hours'.format(t))
    f.close()
