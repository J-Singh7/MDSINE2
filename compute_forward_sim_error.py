'''Compute forward simulation error

Author: David Kaplan
Date: 12/06/20
MDSINE2 version: 4.0.6

Array shapes
------------
N_O : number of taxas
N_T : number of timepoints
N_g : number of Gibb steps

Predicted trajectory
    2-dimensional
        If the array is 2-dim, then we assume the shape is (N_O, N_T)
    3-dimensional
        If the array is 2-dim, then we assume the shape is (N_g, N_O, N_T)
        Error is reported as the mean over the Gibb steps as a single number, 
        the error is not saved for each Gibb step
Ground truth
    The only acceptable input shape is (N_O, N_T)

Error types
-----------
Specify the type of error with `--error` parameter.

Options:
    `relRMSE` - relative Root Mean Square Error (RMSE)
        Convert the input data into relative abudnance and then calculate 
        the RMSE
    `spearman` - Spearman Correlation
        Calculate the spearman correlation for each Taxa's trajectory, then return 
        the mean spearman correlation over all of the taxa
    `RMSE` - Root Mean Square Error
        Compute the straight RMSE

Input Folder
------------
The parameter `--input` is a path of a folder that contains a set of `.npy` files. The 
files are from the output of the script `forward_sim.py`. The files come in pairs:

Full prediction
    {study_name}-{subject_name}-full.py : Predicted
    {study_name}-{subject_name}-full-truth.py : Ground truth
Time lookahead
    {studyname}-{subjname}-start{start}-ndays{ndays}.npy : Predicted
    {studyname}-{subjname}-start{start}-ndays{ndays}-truth.npy : Ground truth

We look for these pairs and automatically load them. If there is not a matching pair
then it will automatically skip over it

Example:
    example_folder/
        healthy-cv2-2-full.npy
        healthy-cv2-2-full-truth.npy
        somethingelse.txt
        healthy-cv3-3-start1.5-ndays7.5.npy
        healthy-cv3-3-start1.5-ndays7.5-truth.npy

    python compute_forward_sim_error.py \\
        --input example_folder \\
        --output example/errors.tsv \\
        --error relRMSE

    This will predict the errors of the pairs
        healthy-cv2-2-full.npy
        healthy-cv2-2-full-truth.npy
    and
        healthy-cv3-3-start1.5-ndays7.5.npy
        healthy-cv3-3-start1.5-ndays7.5-truth.npy
    and will skip over 
        somethingelse.txt

Output
------
The output is a `.tsv` file with the columns:
    study : str
        Name of the study
    subject : str
        Name of the subject
    start : float
        Start day. If a full prediction this is set to zero
    n-days : float
        Number of days looked ahead. If a full prediction this is set to NaN
    error-type : str
        This is the type of error
    error : float
        This is the magnitude of the error
'''

import mdsine2 as md2
import numpy as np
import os
import logging
import argparse
import scipy.stats
import re
import pandas as pd

def _relRMSE_err(pred, truth):
    '''Relative Root Mean Square Error

    Parameters
    ----------
    pred : np.ndarray 2-dim or 3-dim
        Predicted trajectory.
        N_O : number of taxas
        N_T : number of timepoints
        N_g : number of Gibb steps
        2-dimensional
            If the array is 2-dim, then we assume the shape is (N_O, N_T)
        3-dimensional
            If the array is 2-dim, then we assume the shape is (N_g, N_O, N_T)
            Error is reported as the mean over the Gibb steps as a single number, 
            the error is not saved for each Gibb step
    truth : np.ndarray 2-dim
        Ground truth array (N_O, N_T)

    Returns
    -------
    float
    '''
    reltruth = truth/np.sum(truth, axis=0)
    # Convert to 3-dim if it is 2-dim
    if pred.ndim == 2:
        pred = pred[np.newaxis, ...] 
    errors = np.zeros(pred.shape[0])
    for i in range(pred.shape[0]):
        relpred = pred[i]/np.sum(pred[i], axis=0)
        errors[i] = np.sqrt(np.mean(np.square(relpred-reltruth)))
    if logging.root.level == logging.DEBUG:
        logging.debug('relRMSE: {}'.format(errors))
    return np.nanmean(errors)

def _RMSE_err(pred, truth):
    '''Root Mean Square Error

    Parameters
    ----------
    pred : np.ndarray 2-dim or 3-dim
        Predicted trajectory.
        N_O : number of taxas
        N_T : number of timepoints
        N_g : number of Gibb steps
        2-dimensional
            If the array is 2-dim, then we assume the shape is (N_O, N_T)
        3-dimensional
            If the array is 2-dim, then we assume the shape is (N_g, N_O, N_T)
            Error is reported as the mean over the Gibb steps as a single number, 
            the error is not saved for each Gibb step
    truth : np.ndarray 2-dim
        Ground truth array (N_O, N_T)

    Returns
    -------
    float
    '''
    # Convert to 3-dim if it is 2-dim
    if pred.ndim == 2:
        pred = pred[np.newaxis, ...]
    errors = np.zeros(pred.shape[0])
    for i in range(pred.shape[0]):
        errors[i] = np.sqrt(np.mean(np.square(pred[i]-truth)))
    if logging.root.level == logging.DEBUG:
        logging.debug('RMSE: {}'.format(errors))
    return np.nanmean(errors)

def _spearman_err(pred, truth):
    '''Mean Spearman correlation

    Parameters
    ----------
    pred : np.ndarray 2-dim or 3-dim
        Predicted trajectory.
        N_O : number of taxas
        N_T : number of timepoints
        N_g : number of Gibb steps
        2-dimensional
            If the array is 2-dim, then we assume the shape is (N_O, N_T)
        3-dimensional
            If the array is 2-dim, then we assume the shape is (N_g, N_O, N_T)
            Error is reported as the mean over the Gibb steps as a single number, 
            the error is not saved for each Gibb step
    truth : np.ndarray 2-dim
        Ground truth array (N_O, N_T)

    Returns
    -------
    float
    '''
    # Convert to 3-dim if it is 2-dim
    if pred.ndim == 2:
        pred = pred[np.newaxis, ...]
    corr = np.zeros(pred.shape[0])
    for i in range(pred.shape[0]):
        temp = np.zeros(pred.shape[1])
        for aidx in range(pred.shape[1]):
            temp[aidx] = scipy.stats.spearmanr(pred[i,aidx,:], truth[aidx,:])[0]
        corr[i] = np.nanmean(temp)
    if logging.root.level == logging.DEBUG:
        logging.debug('Spearman Correlation: {}'.format(corr))
    return np.nanmean(corr)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--input', '-i', type=str, dest='input_basepath',
        help='Location of the predicted trajectory')
    parser.add_argument('--error', '-e', type=str, dest='errors', nargs='+',
        help='Error Type. You can pass in multiple errors. Options: ' \
            '"relRMSE", "spearman", "RMSE"')
    parser.add_argument('--output', '-o', type=str, dest='output',
        help='Name of the table to save it as')
    parser.add_argument('--sep', '-s', type=str, dest='sep', default='\t',
        help='Separator for the output table')

    args = parser.parse_args()
    md2.LoggingConfig(level=logging.DEBUG)
    input_basepath = args.input_basepath

    # Check if the input path exists
    if not os.path.isdir(input_basepath):
        raise ValueError('Input basepath `{}` does not exist or is not a directory. ' \
            'Check input parameters'.format(input_basepath))

    # Check if the errors are okay
    for error in args.errors:
        if error not in ['relRMSE', 'spearman', 'RMSE']:
            raise ValueError('error ({}) not recognized'.format(error))

    # Iterate over all of the files
    fnames = os.listdir(input_basepath)
    fnames_used = set([])
    data = []
    columns = ['study', 'subject', 'start', 'n-days', 'error-type', 'error']

    re_tla_pred = re.compile(r'^(.*)-(.*)-start(.*)-ndays(.*).npy$')
    re_full_pred = re.compile(r'^(.*)-(.*)-full.npy$')
    
    for error in args.errors:

        if error == 'relRMSE':
            error_func = _relRMSE_err
        elif error == 'RMSE':
            error_func = _RMSE_err
        else:
            error_func = _spearman_err

        logging.info('Error {} with function {}'.format(error, error_func.__name__))

        for fname in fnames:
            if fname in fnames_used:
                continue
            if '-times.npy' in fname:
                # Do not iterate over these (not necessary)
                continue

            if 'full' in  fname:
                # Full time lookahead
                if 'truth' in fname:
                    truth_fname = fname
                    pred_fname = fname.replace('-truth.npy', '.npy')
                else:
                    pred_fname = fname
                    truth_fname = fname.replace('.npy', '-truth.npy')

                # add names to already found
                fnames_used.add(pred_fname)
                fnames_used.add(truth_fname)

                try:
                    studyname, subjectname = re_full_pred.findall(pred_fname)[0]
                except:
                    logging.warning('{} not recognized as a pattern. skipping'.format(fname))
                    continue
                start = np.nan
                ndays = np.nan
            else:
                if 'truth' in fname:
                    truth_fname = fname
                    pred_fname = fname.replace('-truth.npy', '.npy')
                else:
                    truth_fname = fname.replace('.npy', '-truth.npy')
                    pred_fname = fname
                try:
                    studyname, subjectname, start, ndays = re_tla_pred.findall(pred_fname)[0]
                except:
                    logging.warning('{} not recognized as a pattern. skipping'.format(fname))
                    continue

            truth_fname = os.path.abspath(os.path.join(input_basepath, truth_fname))
            pred_fname = os.path.abspath(os.path.join(input_basepath, pred_fname))

            truth = np.load(truth_fname)
            pred = np.load(pred_fname)

            # check the shapes
            if truth.ndim != 2:
                raise ValueError('truth array {} has the shape {}. Not recognized'.format(
                    truth_fname, truth.shape))
            n_taxa_truth = truth.shape[0]
            n_times_truth = truth.shape[1]

            if pred.ndim not in [2,3]:
                raise ValueError('predictin array {} has the shape {}. Not recognized'.format( 
                    pred_fname, pred.shape))
            if pred.ndim == 2:
                n_taxa_pred = pred.shape[0]
                n_times_pred = pred.shape[1]
            else:
                n_taxa_pred = pred.shape[1]
                n_times_pred = pred.shape[2]

            if n_taxa_truth != n_taxa_pred or n_times_pred != n_times_truth:
                raise ValueError('pred shape {} not the same as truth shape {}'.format( 
                    (n_taxa_pred, n_times_pred), (n_taxa_truth, n_times_truth)))
            
            logging.info('Calculating {}-{}-{}-{}'.format(studyname, subjectname, start, ndays))

            err = error_func(pred=pred, truth=truth)
            data.append([studyname, subjectname, start, ndays, error, err])

    # Write the table
    df = pd.DataFrame(data, columns=columns)
    df.to_csv(args.output, sep=args.sep, index=False, header=True)