'''Run leave-one-out cross validation.

Author: David Kaplan
Date: 11/30/20
MDSINE2 version: 4.0.6

This script runs inference for each cross validation for a single fold. Specify the
fold by saying which subject (by name) to leave out
'''
import mdsine2 as md2
import argparse
import logging
import os
import pathlib
import sys

command_fmt = 'python {script} --input {dset} ' \
    '--negbin {negbin} ' \
    '--seed {seed} ' \
    '--burnin {burnin} ' \
    '--n-samples {n_samples} ' \
    '--checkpoint {ckpt} ' \
    '--basepath {basepath} ' \
    '--multiprocessing {mp} ' \
    '--interaction-ind-prior {interaction_prior} ' \
    '--perturbation-ind-prior {perturbation_prior}' 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--dataset', '-d', type=str, dest='dataset',
        help='This is the Gibson dataset we want to do cross validation on')
    parser.add_argument('--cv-basepath', '-o', type=str, dest='output_basepath',
        help='This is the basepath to save the output')
    parser.add_argument('--dset-basepath', '-db', type=str, dest='input_basepath',
        help='This is the basepath to load and save the cv datasets')
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
    parser.add_argument('--basepath', '-b', type=str, dest='basepath',
        help='This is folder to save the output of inference')
    parser.add_argument('--multiprocessing', '-mp', type=int, dest='mp',
        help='If 1, run the inference with multiprocessing. Else run on a single process',
        default=0)
    parser.add_argument('--interaction-ind-prior', '-ip', type=str, dest='interaction_prior',
        help='Prior of the indicator of the interactions')
    parser.add_argument('--perturbation-ind-prior', '-pp', type=str, dest='perturbation_prior',
        help='Prior of the indicator of the perturbations')
    
    args = parser.parse_args()
    md2.config.LoggingConfig(level=logging.INFO)

    input_basepath = args.input_basepath
    os.makedirs(input_basepath, exist_ok=True)
    os.makedirs(args.output_basepath, exist_ok=True)

    logging.info('Loading dataset {}'.format(args.dataset))
    study_master = md2.Study.load(args.dataset)
    subj = study_master[args.leave_out_subj]

    logging.info('Leave out {}'.format(subj.name))
    study = md2.Study.load(args.dataset)
    val_study = study.pop_subject(subj.name)
    study.name = study.name + '-cv{}'.format(subj.name)
    val_study.name = study.name + '-validate'

    # Save the datasets
    study_fname = os.path.join(input_basepath, study.name + '.pkl')
    study.save(study_fname)
    val_study.save(os.path.join(input_basepath, val_study.name + '.pkl'))

    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'step_5_infer_mdsine2.py')
    path = '"' + path + '"'

    logging.info('Run inference')
    command = command_fmt.format(
        script=path,
        dset=study_fname, negbin=args.negbin, seed=args.seed, 
        burnin=args.burnin, n_samples=args.n_samples, ckpt=args.checkpoint,
        basepath=args.output_basepath, mp=args.mp,
        interaction_prior=args.interaction_prior,
        perturbation_prior=args.perturbation_prior)
    logging.info(command)
    os.system(command)