'''Learn the Negative Binomial dispersion parameters that parameterize the noise
of the reads. We learn these parameters with replicate measurements of the reads.

Author: David Kaplan
Date: 11/30/20
MDSINE2 version: 4.0.6
'''
import mdsine2 as md2
import argparse
import logging
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--input', '-i', type=str, dest='input',
        help='This is the dataset to do inference with.')
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
    parser.add_argument('--basepath', '-b', type=str, dest='basepath',
        help='This is folder to save the output of inference')
    args = parser.parse_args()

    study = md2.Study.load(args.input)
    os.makedirs(args.basepath, exist_ok=True)
    basepath = os.path.join(args.basepath, study.name)
    os.makedirs(basepath, exist_ok=True)
    md2.config.LoggingConfig(level=logging.INFO, basepath=basepath)

    # 1) Load the parameters    
    params = md2.config.NegBinConfig(seed=args.seed, burnin=args.burnin, n_samples=args.n_samples,
        ckpt=args.checkpoint, basepath=basepath)
    if args.multiprocessing == 1:
        params.MP_FILTERING = 'full'
    elif args.multiprocessing == 0:
        params.MP_FILTERING = 'debug'
    else:
        raise ValueError('`multiprocessing` ({}) not recognized'.format(args.multiprocessing))

    # 2) Perform inference
    mcmc = md2.negbin.build_graph(params=params, graph_name=study.name, subjset=study)
    mcmc = md2.negbin.run_graph(mcmc, crash_if_error=True)
    mcmc.save()
    study.save(os.path.join(params.MODEL_PATH, md2.config.SUBJSET_FILENAME))

