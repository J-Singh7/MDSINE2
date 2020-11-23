'''Posterior objects used for inference of MDSINE2 model
'''
import logging
import time
import itertools
import psutil
import os

import numpy as np
import numba
import scipy.sparse
import numpy.random as npr
import scipy.stats
import scipy.sparse
import scipy
import math
import random

from .names import STRNAMES, REPRNAMES
from . import pylab as pl

from . import visualization
import matplotlib.pyplot as plt
_LOG_INV_SQRT_2PI = np.log(1/np.sqrt(2*math.pi))

# Helper functions
#-----------------
def _normal_logpdf(value, mean, std):
    '''We use this function if `pylab.random.normal.logpdf` fails to compile,
    which can happen when running jobs on the cluster.
    '''
    return _LOG_INV_SQRT_2PI + (-0.5*((value-mean)/std)**2) - np.log(std)

def negbin_loglikelihood(k,m,dispersion):
    '''Loglikelihood - with parameterization in [1]

    Parameters
    ----------
    k : int
        Observed counts
    m : int
        Mean
    phi : float
        Dispersion

    Returns
    -------
    float
        Negative Binomial Log Likelihood

    References
    ----------
    [1] TE Gibson, GK Gerber. Robust and Scalable Models of Microbiome Dynamics. ICML (2018)
    '''
    r = 1/dispersion
    return math.lgamma(k+r) - math.lgamma(k+1) - math.lgamma(r) \
            + r * (math.log(r) - math.log(r+m)) + k * (math.log(m) - math.log(r+m))

@numba.jit(nopython=True, fastmath=True, cache=False)
def negbin_loglikelihood_MH_condensed(k,m,dispersion):
        '''
        Loglikelihood - with parameterization in [1] - but condensed (do not calculate stuff
        we do not have to)

        Parameters
        ----------
        k : int
            Observed counts
        m : int
            Mean
        phi : float
            Dispersion

        Returns
        -------
        float
            Negative Binomial Log Likelihood

        References
        ----------
        [1] TE Gibson, GK Gerber. Robust and Scalable Models of Microbiome Dynamics. ICML (2018)
        '''
        r = 1/dispersion
        rm = r+m
        return math.lgamma(k+r) - math.lgamma(r) \
            + r * (math.log(r) - math.log(rm)) + k * (math.log(m) - math.log(rm))

def negbin_loglikelihood_MH_condensed_not_fast(k,m,dispersion):
        '''
        Loglikelihood - with parameterization in [1] - but condensed (do not calculate stuff
        we do not have to). We use this function if `negbin_loglikelihood_MH_condensed` fails to
        compile, which can happen when doing jobs on the cluster

        Parameters
        ----------
        k : int
            Observed counts
        m : int
            Mean
        phi : float
            Dispersion

        Returns
        -------
        float
            Negative Binomial Log Likelihood

        References
        ----------
        [1] TE Gibson, GK Gerber. Robust and Scalable Models of Microbiome Dynamics. ICML (2018)
        '''
        r = 1/dispersion
        rm = r+m
        return math.lgamma(k+r) - math.lgamma(r) \
            + r * (math.log(r) - math.log(rm)) + k * (math.log(m) - math.log(rm))

def expected_n_clusters(G):
    '''Calculate the expected number of clusters given the number of ASVs

    Parameters
    ----------
    G : pl.Graph
        Graph object

    Returns
    -------
    int
        Expected number of clusters
    '''
    conc = G[STRNAMES.CONCENTRATION].prior.mean()
    return conc * np.log((G.data.n_asvs + conc) / conc)

def build_prior_covariance(G, cov, order, sparse=True, diag=False, cuda=False):
    '''Build basic prior covariance or precision for the variables
    specified in `order`

    Parameters
    ----------
    G : pylab.graph.Graph
        Graph to get the variables from
    cov : bool
        If True, build the covariance. If False, build the precision
    order : list(str)
        Which parameters to get the priors of
    sparse : bool
        If True, return as a sparse matrix
    diag : bool
        If True, returns the diagonal of the matrix. If this is True, it
        overwhelms the flag `sparse`
    cuda : bool
        If True, returns the array/matrix on the gpu (if there is one). Will not return 
        in sparse form - only dense.

    Returns
    -------
    arr : np.ndarray, scipy.sparse.dia_matrix, torch.DoubleTensor
        Prior covariance or precision matrix in either dense (np.ndarray) or
        sparse (scipy.sparse.dia_matrix) form
    '''
    n_asvs = G.data.n_asvs
    a = []
    for reprname in order:
        if reprname == REPRNAMES.GROWTH_VALUE:
            a.append(np.full(n_asvs, G[REPRNAMES.PRIOR_VAR_GROWTH].value))

        elif reprname == REPRNAMES.SELF_INTERACTION_VALUE:
            a.append(np.full(n_asvs, G[REPRNAMES.PRIOR_VAR_SELF_INTERACTIONS].value))

        elif reprname == REPRNAMES.CLUSTER_INTERACTION_VALUE:
            n_interactions = G[REPRNAMES.CLUSTER_INTERACTION_INDICATOR].num_pos_indicators
            a.append(np.full(n_interactions, G[REPRNAMES.PRIOR_VAR_INTERACTIONS].value))

        elif reprname == REPRNAMES.PERT_VALUE:
            for perturbation in G.perturbations:
                num_on = perturbation.indicator.num_on_clusters()
                a.append(np.full(
                    num_on,
                    perturbation.magnitude.prior.var.value))

        else:
            raise ValueError('reprname ({}) not recognized'.format(reprname))

    if len(a) == 1:
        arr = np.asarray(a[0])
    else:
        arr = np.asarray(list(itertools.chain.from_iterable(a)))
    if not cov:
        arr = 1/arr
    # if cuda:
    #     arr = torch.DoubleTensor(arr).to(_COMPUTE_DEVICE)
    if diag:
        return arr
    # if cuda:
    #     return torch.diag(arr)
    if sparse:
        return scipy.sparse.dia_matrix((arr,[0]), shape=(len(arr),len(arr))).tocsc()
    else:
        return np.diag(arr)

def build_prior_mean(G, order, shape=None, cuda=False):
    '''Builds the prior mean vector for all the variables in `order`.

    Parameters
    ----------
    G : pylab.grapg.Graph
        Graph to index the objects
    order : list
        list of objects to add the priors of. If the variable is the
        cluster interactions or cluster perturbations, then we assume the
        prior mean is a scalar and we set that value for every single value.
    shape : tuple, None
        Shape to cast the array into
    cuda : bool
        If True, returns the array/matrix on the gpu (if there is one)

    Returns
    -------
    np.ndarray, torch.DoubleTensor
    '''
    a = []
    for name in order:
        v = G[name]
        if v.id == REPRNAMES.GROWTH_VALUE:
            a.append(v.prior.mean.value * np.ones(G.data.n_asvs))
        elif v.id == REPRNAMES.SELF_INTERACTION_VALUE:
            a.append(v.prior.mean.value * np.ones(G.data.n_asvs))
        elif v.id == REPRNAMES.CLUSTER_INTERACTION_VALUE:
            a.append(
                np.full(
                    G[REPRNAMES.CLUSTER_INTERACTION_INDICATOR].num_pos_indicators,
                    v.prior.mean.value))
        elif v.id == REPRNAMES.PERT_VALUE:
            for perturbation in G.perturbations:
                a.append(np.full(
                    perturbation.indicator.num_on_clusters(),
                    perturbation.magnitude.prior.mean.value))
        else:
            raise ValueError('`name` ({}) not recognized'.format(name))
    if len(a) == 1:
        a = np.asarray(a[0])
    else:
        a = np.asarray(list(itertools.chain.from_iterable(a)))
    if shape is not None:
        a = a.reshape(*shape)
    # if cuda:
    #     a = torch.DoubleTensor(a).to(_COMPUTE_DEVICE)
    return a

def sample_categorical_log(log_p):
    '''Generate one sample from a categorical distribution with event
    probabilities provided in unnormalized log-space.

    Parameters
    ----------
    log_p : array_like
        logarithms of event probabilities, ***which need not be normalized***

    Returns
    -------
    int
        One sample from the categorical distribution, given as the index of that
        event from log_p.
    '''
    try:
        exp_sample = math.log(random.random())
        events = np.logaddexp.accumulate(np.hstack([[-np.inf], log_p]))
        events -= events[-1]
        return next(x[0]-1 for x in enumerate(events) if x[1] >= exp_sample)
    except:
        logging.critical('CRASHED IN `sample_categorical_log`:\nlog_p{}'.format(
            log_p))
        raise

def log_det(M, var):
    '''Computes pl.math.log_det but also saves the array if it crashes

    Parameters
    ----------
    M : nxn matrix (np.ndarray, scipy.sparse)
        Matrix to calculate the log determinant
    var : pl.variable.Variable subclass
        This is the variable that `log_det` was called from

    Returns
    -------
    np.ndarray
        Log determinant of matrix
    '''
    if scipy.sparse.issparse(M):
        M_ = np.zeros(shape=M.shape)
        M.toarray(out=M_)
        M = M_
    try:
        # if type(M) == torch.Tensor:
        #     return torch.inverse(M)
        # else:
        return pl.math.log_det(M)
    except:
        try:
            sample_iter = var.sample_iter
        except:
            sample_iter = None
        filename = 'crashes/logdet_error_iter{}_var{}pinv_{}.npy'.format(
            sample_iter, var.name, var.G.name)
        logging.critical('\n\n\n\n\n\n\n\nSaved array at "{}" - now crashing\n\n\n'.format(
                filename))
        os.makedirs('crashes/', exist_ok=True)
        np.save(filename, M)
        raise

def pinv(M, var):
    '''Computes np.linalg.pinv but it also saves the array that crashed it if
    it crashes.

    Parameters
    ----------
    M : nxn matrix (np.ndarray, scipy.sparse)
        Matrix to invert
    var : pl.variable.Variable subclass
        This is the variable that `pinv` was called from

    Returns
    -------
    np.ndarray
        Inverse of the matrix
    '''
    if scipy.sparse.issparse(M):
        M_ = np.zeros(shape=M.shape)
        M.toarray(out=M_)
        M = M_
    try:
        # if type(M) == torch.Tensor:
        #     return torch.inverse(M)
        # else:
        try:
            return np.linalg.pinv(M)
        except:
            try:
                return scipy.linalg.pinv(M)
            except:
                return scipy.linalg.inv(M)
    except:
        try:
            sample_iter = var.sample_iter
        except:
            sample_iter = None
        filename = 'crashes/pinv_error_iter{}_var{}pinv_{}.npy'.format(
            sample_iter, var.name, var.G.name)
        logging.critical('\n\n\n\n\n\n\n\nSaved array at "{}" - now crashing\n\n\n'.format(
                filename))
        os.makedirs('crashes/', exist_ok=True)
        np.save(filename, M)
        raise

def _scalar_visualize(obj, path, f, section='posterior'):
    '''Render the traces in the folder `basepath` and write the 
    learned values to the file `f`. This works for scalar variables

    Parameters
    ----------
    obj : mdsine2.Variable
    path : str
        This is the path to write the files to
    f : _io.TextIOWrapper
        File that we are writing the values to
    section : str
        Section of the trace to compute on. Options:
            'posterior' : posterior samples
            'burnin' : burn-in samples
            'entire' : both burn-in and posterior samples

    Returns
    -------
    _io.TextIOWrapper
    '''
    f.write('\n\n###################################\n{}'.format(obj.name))
    f.write('\n###################################\n')
    if not obj.G.inference.is_being_traced(obj):
        f.write('`{}` not learned\n\tValue: {}\n'.format(obj.name, obj.value))
        return f

    summ = pl.summary(obj, section=section)
    for k,v in summ.items():
        f.write('\t{}: {}\n'.format(k,v))

    ax1, _ = visualization.render_trace(var=obj, plt_type='both', 
        section=section, include_burnin=True, log_scale=True, rasterized=True)

    if pl.hasprior(obj):
        
        l,h = ax1.get_xlim()
        try:
            xs = np.arange(l,h,step=(h-l)/1000) 
            ys = []
            for x in xs:
                ys.append(obj.prior.pdf(value=x))
            ax1.plot(xs, ys, label='prior', alpha=0.5, color='red')
            ax1.legend()
        except OverflowError:
            logging.critical('OverflowError while plotting prior')

    fig = plt.gcf()
    fig.suptitle(obj.name)
    plt.savefig(path)
    plt.close()
    return f

# @numba.jit(nopython=True, fastmath=True, cache=True)
def prod_gaussians(means, variances):
    '''Product of Gaussians

    $\mu = [\mu_1, \mu_2, ..., \mu_n]$
    $\var = [\var_1, \var_2, ..., \var_3]$

    Means and variances must be in the same order.

    Parameters
    ----------
    means : np.ndarray
        All of the means
    variances : np.ndarray
        All of the means
    '''
    def _calc_params(mu1, mu2, var1, var2):
        v = var1+var2
        mu = ((var1*mu2) + (var2*mu1))/(v)
        var = (var1*var2)/v
        return mu,var
    mu = means[0]
    var = variances[0]
    for i in range(1,len(means)):
        mu, var = _calc_params(mu1=mu, mu2=means[i], var1=var, var2=variances[i])
    return mu, var

class _Loess(object):
    '''LOESS - Locally Estimated Scatterplot Smoothing
    This module was created by João Paulo Figueira and copied from the 
    repository: https://github.com/joaofig/pyloess.git
    There are a few modifications, mostly to handle edge cases, i.e., what 
    happens when the entire thing is zero
    '''
    def __init__(self, xx, yy, degree=1):
        self.n_xx, self.min_xx, self.max_xx = self.normalize_array(xx)
        self.n_yy, self.min_yy, self.max_yy = self.normalize_array(yy)
        self.degree = degree

        self.input_yy = yy

    @staticmethod
    def _tricubic(x):
        y = np.zeros_like(x)
        idx = (x >= -1) & (x <= 1)
        y[idx] = np.power(1.0 - np.power(np.abs(x[idx]), 3), 3)
        return y

    @staticmethod
    def normalize_array(array):
        min_val = np.min(array)
        max_val = np.max(array)
        return (array - min_val) / (max_val - min_val), min_val, max_val

    @staticmethod
    def get_min_range(distances, window):
        min_idx = np.argmin(distances)
        n = len(distances)
        if min_idx == 0:
            return np.arange(0, window)
        if min_idx == n-1:
            return np.arange(n - window, n)

        min_range = [min_idx]
        while len(min_range) < window:
            i0 = min_range[0]
            i1 = min_range[-1]
            if i0 == 0:
                min_range.append(i1 + 1)
            elif i1 == n-1:
                min_range.insert(0, i0 - 1)
            elif distances[i0-1] < distances[i1+1]:
                min_range.insert(0, i0 - 1)
            else:
                min_range.append(i1 + 1)
        return np.array(min_range)

    @staticmethod
    def get_weights(distances, min_range):
        max_distance = np.max(distances[min_range])
        weights = _Loess._tricubic(distances[min_range] / max_distance)
        return weights

    def normalize_x(self, value):
        return (value - self.min_xx) / (self.max_xx - self.min_xx)

    def denormalize_y(self, value):
        return value * (self.max_yy - self.min_yy) + self.min_yy

    def estimate(self, x, window, use_matrix=False, degree=1):
        n_x = self.normalize_x(x)
        distances = np.abs(self.n_xx - n_x)
        min_range = self.get_min_range(distances, window)
        weights = self.get_weights(distances, min_range)

        if use_matrix or degree > 1:
            wm = np.multiply(np.eye(window), weights)
            xm = np.ones((window, degree + 1))

            xp = np.array([[math.pow(n_x, p)] for p in range(degree + 1)])
            for i in range(1, degree + 1):
                xm[:, i] = np.power(self.n_xx[min_range], i)

            ym = self.n_yy[min_range]
            xmt_wm = np.transpose(xm) @ wm
            beta = np.linalg.pinv(xmt_wm @ xm) @ xmt_wm @ ym
            y = (beta @ xp)[0]
        else:
            xx = self.n_xx[min_range]
            yy = self.n_yy[min_range]
            sum_weight = np.sum(weights)
            sum_weight_x = np.dot(xx, weights)
            sum_weight_y = np.dot(yy, weights)
            sum_weight_x2 = np.dot(np.multiply(xx, xx), weights)
            sum_weight_xy = np.dot(np.multiply(xx, yy), weights)

            mean_x = sum_weight_x / sum_weight
            mean_y = sum_weight_y / sum_weight

            b = (sum_weight_xy - mean_x * mean_y * sum_weight) / \
                (sum_weight_x2 - mean_x * mean_x * sum_weight)
            a = mean_y - b * mean_x
            y = a + b * n_x
        
        ret = self.denormalize_y(y)
        if np.isnan(ret):
            if np.all(self.input_yy == 0):
                return 0
            else:
                raise ValueError('`Returning `np.nan`')
        return ret


# Process Variance
# ----------------
class ProcessVarGlobal(pl.variables.SICS):
    '''Learn a Process variance where we learn th same process variance
    for each ASV. This assumes that the model we're using uses the logscale
    of the data.
    '''
    def __init__(self, prior, **kwargs):
        '''
        Parameters
        ----------
        prior : pl.variables.SICS
            This is the prior of the distribution
        kwargs : dict
            These are the extra parameters for the Variable class
        '''
        kwargs['name'] = STRNAMES.PROCESSVAR
        pl.variables.SICS.__init__(self, dtype=float, **kwargs)
        self.add_prior(prior)
        self.global_variance = True
        self._strr = 'NA'

    def __str__(self):
        return self._strr

    def initialize(self, dof_option, scale_option, value_option, 
        dof=None, scale=None, value=None, variance_scaling=1,
        delay=0):
        '''Initialize the value and hyperparameter.

        Parameters
        ----------
        dof_option : str
            How to initialize the `dof` parameter. Options:
                'manual'
                    Manually specify the dof with the parameter `dof`
                'half'
                    Set the degrees of freedom to the number of data points
                'auto', 'diffuse'
                    Set the degrees of freedom to a sparse number (2.5)
        scale_option : str
            How to initialize the scale of the parameter. Options:
                'manual'
                    Need to also specify `scale` parameter
                'med', 'auto'
                    Set the scale such that mean of the distribution
                    has medium noise (20%)
                'low'
                    Set the scale such that mean of the distribution
                    has low noise (10%)
                'high'
                    Set the scale such that mean of the distribution
                    has high noise (30%)
        value_option : str
            How to initialize the value
                'manual'
                    Set the value with the `value` parameter
                'prior-mean', 'auto'
                    Set the value to the mean of the prior
        variance_scaling : float, None
            How much to inflate the variance
        '''
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        self.delay = delay

        # Set the dof
        if not pl.isstr(dof_option):
            raise TypeError('`dof_option` ({}) must be a str'.format(type(dof_option)))
        if dof_option == 'manual':
            if not pl.isnumeric(dof):
                raise TypeError('`dof` ({}) must be a numeric'.format(type(dof)))
            if dof <= 0:
                raise ValueError('`dof` must be > 0'.format(dof))
            if dof <= 2:
                logging.critical('Process Variance dof ({}) is set unproper'.format(dof))
        elif dof_option == 'half':
            dof = len(self.G.data.lhs)
        elif dof_option in ['auto', 'diffuse']:
            dof = 2.5
        else:
            raise ValueError('`dof_option` ({}) not recognized'.format(dof_option))
        self.prior.dof.override_value(dof)

        # Set the scale
        if not pl.isstr(scale_option):
            raise TypeError('`scale_option` ({}) must be a str'.format(type(scale_option)))
        if scale_option == 'manual':
            if not pl.isnumeric(scale):
                raise TypeError('`scale` ({}) must be a numeric'.format(type(scale)))
            if scale <= 0:
                raise ValueError('`scale` ({}) must be > 0'.format(scale))
        elif scale_option in ['auto', 'med']:
            scale = (0.2 ** 2) * (self.prior.dof.value - 2) / self.prior.dof.value
        elif scale_option ==  'low':
            scale = (0.1 ** 2) * (self.prior.dof.value - 2) / self.prior.dof.value
        elif scale_option ==  'high':
            scale = (0.3 ** 2) * (self.prior.dof.value - 2) / self.prior.dof.value
        else:
            raise ValueError('`scale_option` ({}) not recognized'.format(scale_option))
        self.prior.scale.override_value(scale)

        # Set the value
        if not pl.isstr(value_option):
            raise TypeError('`value_option` ({}) must be a str'.format(type(value_option)))
        if value_option == 'manual':
            if not pl.isnumeric(dof):
                raise TypeError('`value` ({}) must be a numeric'.format(type(value)))
            if value <= 0:
                raise ValueError('`value` ({}) must be > 0'.format(value))
        elif value_option in ['prior-mean', 'auto']:
            value = self.prior.mean()
        else:
            raise ValueError('`value_option` ({}) not recognized'.format(value_option))
        self.value = value
        
        self.rebuild_diag()
        self._there_are_perturbations = self.G.perturbations is not None

    def build_matrix(self, cov, sparse=True):
        '''Builds the process variance as a covariance or precision
        matrix.

        Parameters
        ----------
        cov : bool
            If True, make the covariance matrix
        sparse : bool
            If True, return the matrix as a sparse matrix

        Returns
        -------
        np.ndarray or scipy.sparse
            Either sparse or dense covariance/precision matrix
        '''
        a = self.diag
        if not cov:
            a = 1/a
        if sparse:
            return scipy.sparse.dia_matrix((a,[0]), shape=(len(a),len(a))).tocsc()
        else:
            return np.diag(a)

    def rebuild_diag(self):
        '''Builds up the process variance diagonal that we use to make the matrix
        '''
        a = self.value / self.G.data.dt_vec
        if self.G.data.zero_inflation_transition_policy is not None:
            a = a[self.G.data.rows_to_include_zero_inflation]
        
        self.diag = a
        self.prec = 1/a

    def update(self):
        '''Update the process variance

        y = (log(x_{k+1}) - log(x_k))/dt
        Xb = a_1 (1 + \gamma) + A x

        % These are our dynamics with the process variance
        y ~ Normal(Xb , \sigma^2_w / dt)

        % Subtract the mean
        y - Xb ~ Normal( 0, \sigma^2_w / dt)
        
        % Substitute
        z = y - Xb 

        z ~ Normal (0, \sigma^2_w / dt)
        z * \sqrt{dt} ~ Normal (0, \sigma^2_w)

        This is now in a form we can use to calculate the posterior
        '''
        if self._there_are_perturbations:
            lhs = [
                REPRNAMES.GROWTH_VALUE, 
                REPRNAMES.SELF_INTERACTION_VALUE,
                REPRNAMES.CLUSTER_INTERACTION_VALUE]
        else:
            lhs = [
                REPRNAMES.GROWTH_VALUE, 
                REPRNAMES.SELF_INTERACTION_VALUE,
                REPRNAMES.CLUSTER_INTERACTION_VALUE]
        
        # This is the residual z = y - Xb
        z = self.G.data.construct_lhs(lhs, 
            kwargs_dict={REPRNAMES.GROWTH_VALUE:{
                'with_perturbations': self._there_are_perturbations}})
        z = np.asarray(z).ravel()
        if self.G.data.zero_inflation_transition_policy is not None:
            z = z * self.G.data.sqrt_dt_vec[self.G.data.rows_to_include_zero_inflation]
        else:
            z = z * self.G.data.sqrt_dt_vec
        residual = np.sum(np.square(z))

        self.dof.value = self.prior.dof.value + len(z)
        self.scale.value = ((self.prior.scale.value * self.prior.dof.value) + \
           residual)/self.dof.value
        
        # shape = 2 + (len(residual)/2)
        # scale = 0.00001 + np.sum(np.square(residual))/2
        # self.value = pl.random.invgamma.sample(shape=shape, scale=scale)
        self.sample()
        self.rebuild_diag()

        self._strr = '{}, empirical_variance: {:.5f}'.format(self.value, 
            residual/len(z))

    def visualize(self, path, f, section='posterior'):
        return _scalar_visualize(self, path=path, f=f, section=section)

# Clustering
# ----------
class Concentration(pl.variables.Gamma):
    '''Defines the posterior for the concentration parameter that is used
    in learning the cluster assignments.
    The posterior is implemented as it is describes in 'Bayesian Inference
    for Density Estimation' by M. D. Escobar and M. West, 1995.
    '''
    def __init__(self, prior, value=None, n_iter=None, **kwargs):
        '''Parameters

        value (float, int)
            - Initial value of the concentration
            - Default value is the mean of the prior
        '''
        kwargs['name'] = STRNAMES.CONCENTRATION
        # initialize shape and scale as the same as the priors
        # we will be updating this later
        pl.variables.Gamma.__init__(self, shape=prior.shape.value, scale=prior.scale.value,
            dtype=float, **kwargs)
        self.add_prior(prior)

    def initialize(self, value_option, hyperparam_option, n_iter=None, value=None,
        shape=None, scale=None, delay=0):
        '''Initialize the hyperparameters of the beta prior

        Parameters
        ----------
        value_option (str)
            - Options to initialize the value
            - 'manual'
                - Set the value manually, `value` must also be specified
            - 'auto', 'prior-mean'
                - Set to the mean of the prior
        hyperparam_option (str)
            - Options ot initialize the hyperparameters
            - Options
                - 'manual'
                    - Set the values manually. `shape` and `scale` must also be specified
                - 'auto', 'diffuse'
                    - shape = 1e-5, scale= 1e5
        shape, scale (int, float)
            - User specified values
            - Only necessary if `hyperparam_option` == 'manual'
        '''
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        self.delay = delay

        if hyperparam_option == 'manual':
            self.prior.shape.override_value(shape)
            self.prior.scale.override_value(scale)
            self.n_iter = n_iter

        elif hyperparam_option in ['diffuse', 'auto']:
            self.prior.shape.override_value(1e-5)
            self.prior.scale.override_value(1e5)
            self.n_iter = 20
        else:
            raise ValueError('hyperparam_option `{}` not recognized'.format(hyperparam_option))

        if value_option == 'manual':
            self.value = value
        elif value_option in ['auto', 'prior-mean']:
            self.value = self.prior.mean()
        else:
            raise ValueError('value_option `{}` not recognized'.format(value_option))

        self.shape.value = self.prior.shape.value
        self.scale.value = self.prior.scale.value
        logging.info('Cluster Concentration initialization results:\n' \
            '\tprior shape: {}\n\tprior scale: {}\n\tvalue: {}'.format(
                self.prior.shape.value, self.prior.scale.value, self.value))

    def update(self):
        '''Sample the posterior of the concentration parameter
        '''
        if self.sample_iter < self.delay:
            return

        clustering = self.G[REPRNAMES.CLUSTER_INTERACTION_VALUE].clustering
        k = len(clustering)
        n = self.G.data.n_asvs
        for i in range(self.n_iter):
            #first sample eta from a beta distribution
            eta = npr.beta(self.value+1,n)
            #sample alpha from a mixture of gammas
            pi_eta = [0.0, n]
            pi_eta[0] = (self.prior.shape.value+k-1.)/(1/(self.prior.scale.value)-np.log(eta))
            self.scale.value =  1/(1/self.prior.scale.value - np.log(eta))
            self.shape.value = self.prior.shape.value + k
            if np.random.choice([0,1], p=pi_eta/np.sum(pi_eta)) != 0:
                self.shape.value -= 1
            self.sample()

    def visualize(self, path, f, section='posterior'):
        '''Render the traces in the folder `basepath` and write the 
        learned values to the file `f`.

        Parameters
        ----------
        path : str
            This is the path to write the files to
        f : _io.TextIOWrapper
            File that we are writing the values to
        section : str
            Section of the trace to compute on. Options:
                'posterior' : posterior samples
                'burnin' : burn-in samples
                'entire' : both burn-in and posterior samples

        Returns
        -------
        _io.TextIOWrapper
        '''
        return _scalar_visualize(self, path=path, f=f, section=section)


class ClusterAssignments(pl.graph.Node):
    '''This is the posterior of the cluster assignments for each ASV.

    To calculate the loglikelihood of an ASV being in a cluster, we have to
    marginalize out the cluster that the ASV belongs to - this means marginalizing
    out the interactions going in and out of the cluster in question, along with all
    the perturbations associated with it.
    '''
    def __init__(self, clustering, concentration, m=1, mp=None, **kwargs):
        '''Parameters

        clustering (pylab.cluster.Clustering)
            - Defines the clusters
        concentration (pylab.Variable or subclass of pylab.Variable)
            - Defines the concentration parameter for the base distribution
        m (int, Optional)
            - Number of auxiliary variables defined in the model
            - Default is 1
        mp : str, None
            This is the type of multiprocessing it is going to be. Options:
                None
                    No multiprocessing
                'full-#'
                    Send out to the different processors, where '#' is the number of
                    processors to make
                'debug'
                    Send out to the different classes but stay on a single core. This
                    is necessary for benchmarking and easier debugging.
        '''
        self.clustering = clustering
        self.concentration = concentration
        self.m = m
        self.mp = mp
        self._strtime = -1

        if self.mp is not None:
            self.update = self.update_mp
        else:
            self.update = self.update_slow_fast

        kwargs['name'] = STRNAMES.CLUSTERING
        pl.graph.Node.__init__(self, **kwargs)

    def __str__(self):
        return str(self.clustering) + '\nTotal time: {}'.format(self._strtime)

    def __getstate__(self):
        state = self.__dict__.copy()
        state.pop('pool', None)
        state.pop('actors', None)
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.pool = []
        self.actors = None

    @property
    def value(self):
        '''This is so that the cluster assignments get printed in inference.MCMC
        '''
        return self.clustering

    @property
    def sample_iter(self):
        '''Sample iteration
        '''
        return self.clustering.n_clusters.sample_iter

    def initialize(self, value_option, hyperparam_option=None, value=None, n_clusters=None,
        delay=0, run_every_n_iterations=1):
        '''Initialize the cluster assingments - there are no hyperparamters to
        initialize because the concentration is initialized somewhere else

        Note - if `n_clusters` is not specified and the cluster initialization
        method requires it - it will be set to the expected number of clusters
        which = log(n_asvs)/log(2)

        Parameters
        ----------
        value_option : str
            The different methods to initialize the clusters
            Options
                'manual'
                    Manually set the cluster assignments
                'no-clusters'
                    Every ASV in their own cluster
                'random'
                    Every ASV is randomly assigned to the number of clusters. `n_clusters` required
                'taxonomy'
                    Cluster ASVs based on their taxonomic similarity. `n_clusters` required
                'sequence'
                    Cluster ASVs based on their sequence similarity. `n_clusters` required
                'phylogeny'
                    Cluster ASVs based on their phylogenetic similarity. `n_clusters` required
                'spearman', 'auto'
                    Creates a distance matrix based on the spearman rank similarity
                    between two trajectories. We use the raw data. `n_clusters` required
                'fixed-topology'
                    Sets the clustering assignment to the most likely clustering configuration
                    specified in the graph at the location `value` (`value` is a str).
                    We take the mean coclusterings and do agglomerative clustering on that matrix
                    with the `mode` number of clusters.
        hyperparam_option : None
            Not used in this function - only here for API consistency
        value : list of list
            Cluster assingments for each of the ASVs
            Only necessary if `value_option` == 'manual'
        n_clusters : int, str
            Necessary if `value_option` is not 'manual' or 'no-clusters'
            If str, options:
                'expected', 'auto': log_2(n_asvs)
        run_every_n_iterations : int
            Only run the update every `run_every_n_iterations` iterations
        '''
        from sklearn.cluster import AgglomerativeClustering
        from .util import generate_cluster_assignments_posthoc
        asvs = self.G.data.asvs

        self.run_every_n_iterations = run_every_n_iterations
        self.delay = delay

        if value_option not in ['manual', 'no-clusters', 'fixed-topology']:
            if pl.isstr(n_clusters):
                if n_clusters in ['expected', 'auto']:
                    n_clusters = expected_n_clusters(self.G)
                else:
                    raise ValueError('`n_clusters` ({}) not recognized'.format(n_clusters))
            if not pl.isint(n_clusters):
                raise TypeError('`n_clusters` ({}) must be a str or an int'.format(type(n_clusters)))
            if n_clusters <= 0:
                raise ValueError('`n_clusters` ({}) must be > 0'.format(n_clusters))
            if n_clusters > self.G.data.n_asvs:
                raise ValueError('`n_clusters` ({}) must be <= than the number of ASVs ({})'.format(
                    n_clusters, self.G.data.n_asvs))

        if value_option == 'manual':
            # Check that all of the ASVs are in the init and that it is in the right structure
            if not pl.isarray(value):
                raise ValueError('if `value_option` is "manual", value ({}) must ' \
                    'be of type array'.format(value.__class__))
            clusters = list(value)

            idxs_to_delete = []
            for idx, cluster in enumerate(clusters):
                if not pl.isarray(cluster):
                    raise ValueError('cluster at index `{}` ({}) is not an array'.format(
                        idx, cluster))
                cluster = list(cluster)
                if len(cluster) == 0:
                    logging.warning('Cluster index {} has 0 elements, deleting'.format(
                        idx))
                    idxs_to_delete.append(idx)
            if len(idxs_to_delete) > 0:
                clusters = np.delete(clusters, idxs_to_delete).tolist()

            all_oidxs = set()
            for cluster in clusters:
                for oidx in cluster:
                    if not pl.isint(oidx):
                        raise ValueError('`oidx` ({}) must be an int'.format(oidx.__class__))
                    if oidx >= len(asvs):
                        raise ValueError('oidx `{}` not in our ASVSet'.format(oidx))
                    all_oidxs.add(oidx)

            for oidx in range(len(asvs)):
                if oidx not in all_oidxs:
                    raise ValueError('oidx `{}` in ASVSet not in `value` ({})'.format(
                        oidx, value))
            # Now everything is checked and valid

        elif value_option == 'fixed-topology':
            logging.info('Fixed topology initialization')
            if not pl.isstr(value):
                raise TypeError('`value` ({}) must be a str'.format(value))

            CHAIN2 = pl.inference.BaseMCMC.load(value)
            CLUSTERING2 = CHAIN2.graph[STRNAMES.CLUSTERING_OBJ]
            ASVS2 = CHAIN2.graph.data.asvs
            asvs_curr = self.G.data.asvs
            for asv in ASVS2:
                if asv.name not in asvs_curr:
                    raise ValueError('Cannot perform fixed topology because the ASV {} in ' \
                        'the passed in clustering is not in this clustering: {}'.format(
                            asv.name, asvs_curr.names.order))
            for asv in asvs_curr:
                if asv.name not in ASVS2:
                    raise ValueError('Cannot perform fixed topology because the ASV {} in ' \
                        'the current clustering is not in the passed in clustering: {}'.format(
                            asv.name, ASVS2.names.order))

            # Get the most likely cluster configuration and set as the value for the passed in cluster
            ret = generate_cluster_assignments_posthoc(CLUSTERING2, n_clusters='mode', set_as_value=False)
            CLUSTERING2.from_array(ret)
            logging.info('Clustering set to:\n{}'.format(str(CLUSTERING2)))

            # Set the passed in cluster assignment as the current cluster assignment
            # Need to be careful because the indices of the ASVs might not line up
            clusters = []
            for cluster in CLUSTERING2:
                anames = [asvs_curr[ASVS2.names.order[aidx]].name for aidx in cluster.members]
                aidxs = [asvs_curr[aname].idx for aname in anames]
                clusters.append(aidxs)

        elif value_option == 'no-clusters':
            clusters = []
            for oidx in range(len(asvs)):
                clusters.append([oidx])

        elif value_option == 'random':
            clusters = {}
            for oidx in range(len(asvs)):
                idx = npr.choice(n_clusters)
                if idx in clusters:
                    clusters[idx].append(oidx)
                else:
                    clusters[idx] = [oidx]
            c = []
            for cid in clusters.keys():
                c.append(clusters[cid])
            clusters = c

        elif value_option == 'taxonomy':
            # Create an affinity matrix, we can precompute the self-similarity to 1
            M = np.diag(np.ones(len(asvs), dtype=float))
            for i, oid1 in enumerate(asvs.ids.order):
                for j, oid2 in enumerate(asvs.ids.order):
                    if i == j:
                        continue
                    M[i,j] = asvs.taxonomic_similarity(oid1=oid1, oid2=oid2)

            c = AgglomerativeClustering(
                n_clusters=n_clusters,
                affinity='precomputed',
                linkage='complete')
            assignments = c.fit_predict(1-M)

            # Convert assignments into clusters
            clusters = {}
            for oidx,cidx in enumerate(assignments):
                if cidx not in clusters:
                    clusters[cidx] = []
                clusters[cidx].append(oidx)
            clusters = [val for val in clusters.values()]

        elif value_option == 'sequence':
            import diversity

            logging.info('Making affinity matrix from sequences')
            evenness = np.diag(np.ones(len(self.G.data.asvs), dtype=float))

            for i in range(len(self.G.data.asvs)):
                for j in range(len(self.G.data.asvs)):
                    if j <= i:
                        continue
                    # Subtract because we want to make a similarity matrix
                    dist = 1-diversity.beta.hamming(
                        list(self.G.data.asvs[i].sequence),
                        list(self.G.data.asvs[j].sequence))
                    evenness[i,j] = dist
                    evenness[j,i] = dist

            c = AgglomerativeClustering(
                n_clusters=n_clusters,
                affinity='precomputed',
                linkage='average')
            assignments = c.fit_predict(evenness)
            clusters = {}
            for oidx,cidx in enumerate(assignments):
                if cidx not in clusters:
                    clusters[cidx] = []
                clusters[cidx].append(oidx)
            clusters = [val for val in clusters.values()]

        elif value_option == 'spearman':
            # Use spearman correlation to create a distance matrix
            # Use agglomerative clustering to make the clusters based
            # on distance matrix (distance = 1 - pearson(x,y))
            dm = np.zeros(shape=(len(asvs), len(asvs)))
            data = []
            for ridx in range(self.G.data.n_replicates):
                data.append(self.G.data.abs_data[ridx])
            data = np.hstack(data)
            for i in range(len(asvs)):
                for j in range(i+1):
                    distance = (1 - scipy.stats.spearmanr(data[i, :], data[j, :])[0])/2
                    dm[i,j] = distance
                    dm[j,i] = distance

            c = AgglomerativeClustering(
                n_clusters=n_clusters,
                affinity='precomputed',
                linkage='complete')
            assignments = c.fit_predict(dm)

            # convert into clusters
            clusters = {}
            for oidx, cidx in enumerate(assignments):
                if cidx not in clusters:
                    clusters[cidx] = []
                clusters[cidx].append(oidx)
            clusters = [val for val in clusters.values()]

        elif value_option == 'phylogeny':
            raise NotImplementedError('`phylogeny` not implemented yet')

        else:
            raise ValueError('`value_option` "{}" not recognized'.format(value_option))

        # Move all the ASVs into their assigned clusters
        for cluster in clusters:
            cid = None
            for oidx in cluster:
                if cid is None:
                    # make new cluster
                    cid = self.clustering.make_new_cluster_with(idx=oidx)
                else:
                    self.clustering.move_item(idx=oidx, cid=cid)
        logging.info('Cluster Assingments initialization results:\n{}'.format(
            str(self.clustering)))
        self._there_are_perturbations = self.G.perturbations is not None

        # Initialize the multiprocessors if necessary
        if self.mp is not None:
            if not pl.isstr(self.mp):
                raise TypeError('`mp` ({}) must be a str'.format(type(self.mp)))
            if 'full' in self.mp:
                n_cpus = self.mp.split('-')[1]
                if n_cpus == 'auto':
                    self.n_cpus = psutil.cpu_count(logical=False)
                else:
                    try:
                        self.n_cpus = int(n_cpus)
                    except:
                        raise ValueError('`mp` ({}) not recognized'.format(self.mp))
                self.pool = pl.multiprocessing.PersistentPool(ptype='dasw', G=self.G)
            elif self.mp == 'debug':
                self.pool = None
            else:
                raise ValueError('`mp` ({}) not recognized'.format(self.mp))
        else:
            self.pool = None

        self.ndts_bias = []
        self.n_asvs = len(self.G.data.asvs)
        self.n_replicates = self.G.data.n_replicates
        self.n_dts_for_replicate = self.G.data.n_dts_for_replicate
        self.total_dts = np.sum(self.n_dts_for_replicate)
        for ridx in range(self.G.data.n_replicates):
            self.ndts_bias.append(
                np.arange(0, self.G.data.n_dts_for_replicate[ridx] * self.n_asvs, self.n_asvs))
        self.replicate_bias = np.zeros(self.n_replicates, dtype=int)
        for ridx in range(1, self.n_replicates):
            self.replicate_bias[ridx] = self.replicate_bias[ridx-1] + \
                self.n_asvs * self.n_dts_for_replicate[ridx - 1]

    def visualize(self, basepath, f, section='posterior', 
        asv_formatter='%(paperformat)s', yticklabels='%(paperformat)s %(index)s', 
        xticklabels='%(index)s'):
        '''Render the traces in the folder `basepath` and write the 
        learned values to the file `f`.

        Parameters
        ----------
        basepath : str
            This is the loction to write the files to
        f : _io.TextIOWrapper
            File that we are writing the values to
        section : str
            Section of the trace to compute on. Options:
                'posterior' : posterior samples
                'burnin' : burn-in samples
                'entire' : both burn-in and posterior samples

        Returns
        -------
        _io.TextIOWrapper
        '''
        from .util import generate_cluster_assignments_posthoc
        asvs = self.G.data.asvs
        f.write('\n\n###################################\n')
        f.write(self.name)
        f.write('\n###################################\n')
        if not self.G.inference.is_in_inference_order(self):
            f.write('`{}` not learned. These were the fixed cluster assignments\n'.format(self.name))
            for cidx, cluster in enumerate(self.clustering):
                f.write('Cluster {}:\n'.format(cidx+1))
                for aidx in cluster.members:
                    label = pl.asvname_formatter(format=asv_formatter, asv=asvs[aidx], asvs=asvs)
                    f.write('\t- {}\n'.format(label))

            return f

        # Coclusters
        cocluster_trace = self.clustering.coclusters.get_trace_from_disk(section=section)
        coclusters = pl.variables.summary(cocluster_trace, section=section)['mean']
        for i in range(coclusters.shape[0]):
            coclusters[i,i] = np.nan

        visualization.render_cocluster_proportions(
            coclusters=coclusters, asvs=self.G.data.asvs, clustering=self.clustering,
            yticklabels=yticklabels, include_tick_marks=False, xticklabels=xticklabels,
            title='Cluster Assignments')
        fig = plt.gcf()
        fig.tight_layout()
        plt.savefig(os.path.join(basepath, 'coclusters.pdf'))
        plt.close()

        # N clusters
        visualization.render_trace(var=self.clustering.n_clusters, plt_type='both', 
            section=section, include_burnin=True, rasterized=True)
        fig = plt.gcf()
        fig.suptitle('Number of Clusters')
        plt.savefig(os.path.join(basepath, 'n_clusters.pdf'))
        plt.close()

        ca = generate_cluster_assignments_posthoc(clustering=self.clustering, n_clusters='mode', 
            section=section)
        cluster_assignments = {}
        for idx, assignment in enumerate(ca):
            if assignment in cluster_assignments:
                cluster_assignments[assignment].append(idx)
            else:
                cluster_assignments[assignment] = [idx]

        f.write('Mode number of clusters: {}\n'.format(len(self.clustering)))
        for idx,lst in enumerate(cluster_assignments.values()):
            f.write('Cluster {} - Size {}\n'.format(idx+1, len(lst)))
            for oidx in lst:
                # Get rid of index because that does not really make sense here
                label = pl.asvname_formatter(format=asv_formatter, asv=asvs[oidx], asvs=asvs)
                f.write('\t- {}\n'.format(label))
        
        return f

    def set_trace(self):
        self.clustering.set_trace()

    def remove_local_trace(self):
        '''Delete the local trace
        '''
        self.clustering.trace = None

    def add_trace(self):
        self.clustering.add_trace()

    def kill(self):
        if pl.ispersistentpool(self.pool):
            # For pylab multiprocessing, explicitly kill them
            self.pool.kill()
        return

    # Update super safe - meant to be used during debugging
    # =====================================================
    def update_slow(self):
        ''' This is updating the new cluster. Depending on the iteration you do
        either split-merge Metropolis-Hasting update or a regular Gibbs update. To
        get highest mixing we alternate between each.
        '''
        if self.clustering.n_clusters.sample_iter < self.delay:
            return

        if self.clustering.n_clusters.sample_iter % self.run_every_n_iterations != 0:
           return

        # print('in clustering')
        start_time = time.time()
        oidxs = npr.permutation(np.arange(len(self.G.data.asvs)))

        for oidx in oidxs:
            self.gibbs_update_single_asv_slow(oidx=oidx)
        self._strtime = time.time() - start_time

    def gibbs_update_single_asv_slow(self, oidx):
        '''The update function is based off of Algorithm 8 in 'Markov Chain
        Sampling Methods for Dirichlet Process Mixture Models' by Radford M.
        Neal, 2000.

        Calculate the marginal likelihood of the asv in every cluster
        and a new cluster then sample from `self.sample_categorical_log`
        to get the cluster assignment.

        Parameters
        ----------
        oidx : int
            ASV index that we are updating the cluster assignment of
        '''
        curr_cluster = self.clustering.idx2cid[oidx]
        concentration = self.concentration.value

        # start as a dictionary then send values to `sample_categorical_log`
        LOG_P = []
        LOG_KEYS = []

        # Calculate current cluster
        # =========================
        # If the element is already in its own cluster, use the new cluster case
        if self.clustering.clusters[curr_cluster].size == 1:
            a = np.log(concentration/self.m)
        else:
            a = np.log(self.clustering.clusters[curr_cluster].size - 1)
        LOG_P.append(a + self.calculate_marginal_loglikelihood_slow()['ret'])
        LOG_KEYS.append(curr_cluster)

        # Calculate going to every other cluster
        # ======================================
        for cid in self.clustering.order:
            if curr_cluster == cid:
                continue

            # Move ASV and recompute the matrices
            self.clustering.move_item(idx=oidx,cid=cid)
            self.G[REPRNAMES.CLUSTER_INTERACTION_INDICATOR].update_cnt_indicators()
            self.G.data.design_matrices[REPRNAMES.CLUSTER_INTERACTION_VALUE].M.build()
            if self._there_are_perturbations:
                self.G.data.design_matrices[REPRNAMES.PERT_VALUE].M.build()

            LOG_P.append(np.log(self.clustering.clusters[cid].size - 1) + \
                self.calculate_marginal_loglikelihood_slow()['ret'])
            LOG_KEYS.append(cid)


        # Calculate new cluster
        # =====================
        cid=self.clustering.make_new_cluster_with(idx=oidx)
        self.G[REPRNAMES.CLUSTER_INTERACTION_INDICATOR].update_cnt_indicators()
        self.G.data.design_matrices[REPRNAMES.CLUSTER_INTERACTION_VALUE].M.build()
        if self._there_are_perturbations:
            self.G.data.design_matrices[REPRNAMES.PERT_VALUE].M.build()
        
        LOG_KEYS.append(cid)
        LOG_P.append(np.log(concentration/self.m) + \
            self.calculate_marginal_loglikelihood_slow()['ret'])

        # Sample the assignment
        # =====================
        idx = sample_categorical_log(LOG_P)
        assigned_cid = LOG_KEYS[idx]
        curr_clus = self.clustering.idx2cid[oidx]

        if assigned_cid != curr_clus:
            self.clustering.move_item(idx=oidx,cid=assigned_cid)
            self.G[REPRNAMES.CLUSTER_INTERACTION_INDICATOR].update_cnt_indicators()

            # Change the mixing matrix for the interactions and (potentially) perturbations
            self.G.data.design_matrices[REPRNAMES.CLUSTER_INTERACTION_VALUE].M.build()
            if self._there_are_perturbations:
                self.G.data.design_matrices[REPRNAMES.PERT_VALUE].M.build()

    def calculate_marginal_loglikelihood_slow(self):
        '''Marginalizes out the interactions and the perturbations
        '''
        # Build the parameters
        # ====================
        self.G[REPRNAMES.CLUSTER_INTERACTION_INDICATOR].update_cnt_indicators()
        lhs = [REPRNAMES.GROWTH_VALUE, REPRNAMES.SELF_INTERACTION_VALUE]
        if self._there_are_perturbations:
            rhs = [REPRNAMES.PERT_VALUE, REPRNAMES.CLUSTER_INTERACTION_VALUE]
        else:
            rhs = [REPRNAMES.CLUSTER_INTERACTION_VALUE]

        # reconstruct the X matrices
        for v in rhs:
            self.G.data.design_matrices[v].M.build()
        
        y = self.G.data.construct_lhs(lhs, 
            kwargs_dict={REPRNAMES.GROWTH_VALUE:{'with_perturbations': False}})
        X = self.G.data.construct_rhs(keys=rhs, toarray=True)
        process_prec = self.G[REPRNAMES.PROCESSVAR].build_matrix(cov=False, sparse=False)
        prior_prec = build_prior_covariance(G=self.G, cov=False, order=rhs, sparse=False)
        prior_var = build_prior_covariance(G=self.G, cov=True, order=rhs, sparse=False)
        prior_mean = build_prior_mean(G=self.G, order=rhs, shape=(-1,1))

        # If nothing is on, return 0
        if X.shape[1] == 0:
            return {
                'a': 0,
                'beta_prec': 0,
                'process_prec': prior_prec,
                'ret': 0,
                'beta_logdet': 0,
                'priorvar_logdet': 0,
                'bEb': 0,
                'bEbprior': 0}

        # Calculate the marginalization
        # =============================
        beta_prec = X.T @ process_prec @ X + prior_prec
        beta_cov = pinv(beta_prec, self)
        beta_mean = beta_cov @ ( X.T @ process_prec @ y + prior_prec @ prior_mean )
        beta_mean = np.asarray(beta_mean).reshape(-1,1)

        try:
            beta_logdet = log_det(beta_cov, self)
        except:
            logging.critical('Crashed in log_det')
            logging.critical('beta_cov:\n{}'.format(beta_cov))
            logging.critical('prior_prec\n{}'.format(prior_prec))
            raise
        priorvar_logdet = log_det(prior_var, self)
        ll2 = 0.5 * (beta_logdet - priorvar_logdet)

        bEbprior = np.asarray(prior_mean.T @ prior_prec @ prior_mean)[0,0]
        bEb = np.asarray(beta_mean.T @ beta_prec @ beta_mean)[0,0]
        ll3 = 0.5 * (bEb  - bEbprior)

        # print('prior_prec truth:\n', prior_prec)

        return {
            'a': X.T @ process_prec,
            'beta_prec': beta_prec,
            'process_prec': prior_prec,
            'ret': ll2+ll3,
            'beta_logdet': beta_logdet,
            'priorvar_logdet': priorvar_logdet,
            'bEb': bEb,
            'bEbprior': bEbprior}

    # Update regular - meant to be used during inference
    # ==================================================
    # @profile
    def update_slow_fast(self):
        '''Much faster than `update_slow`
        '''

        if self.clustering.n_clusters.sample_iter < self.delay:
            return

        if self.clustering.n_clusters.sample_iter % self.run_every_n_iterations != 0:
           return

        start_time = time.time()

        self.process_prec = self.G[REPRNAMES.PROCESSVAR].prec.ravel() #.build_matrix(cov=False, sparse=False)
        self.process_prec_matrix = self.G[REPRNAMES.PROCESSVAR].build_matrix(sparse=True, cov=False)
        lhs = [REPRNAMES.GROWTH_VALUE, REPRNAMES.SELF_INTERACTION_VALUE]
        self.y = self.G.data.construct_lhs(lhs, 
            kwargs_dict={REPRNAMES.GROWTH_VALUE:{'with_perturbations': False}})

        oidxs = npr.permutation(np.arange(len(self.G.data.asvs)))
        iii = 0
        for oidx in oidxs:
            logging.info('{}/{}: {}'.format(iii, len(oidxs), oidx))
            self.gibbs_update_single_asv_slow_fast(oidx=oidx)
            iii += 1
        self._strtime = time.time() - start_time
    
    # @profile
    def gibbs_update_single_asv_slow_fast(self, oidx):
        '''The update function is based off of Algorithm 8 in 'Markov Chain
        Sampling Methods for Dirichlet Process Mixture Models' by Radford M.
        Neal, 2000.

        Calculate the marginal likelihood of the asv in every cluster
        and a new cluster then sample from `self.sample_categorical_log`
        to get the cluster assignment.

        Parameters
        ----------
        oidx : int
            ASV index that we are updating the cluster assignment of
        '''
        curr_cluster = self.clustering.idx2cid[oidx]
        concentration = self.concentration.value

        # start as a dictionary then send values to `sample_categorical_log`
        LOG_P = []
        LOG_KEYS = []

        # Calculate current cluster
        # =========================
        # If the element is already in its own cluster, use the new cluster case
        if self.clustering.clusters[curr_cluster].size == 1:
            a = np.log(concentration/self.m)
        else:
            a = np.log(self.clustering.clusters[curr_cluster].size - 1)
        LOG_P.append(a + self.calculate_marginal_loglikelihood_slow_fast_sparse())
        LOG_KEYS.append(curr_cluster)

        # Calculate going to every other cluster
        # ======================================
        for cid in self.clustering.order:
            if curr_cluster == cid:
                continue

            # Move ASV and recompute the matrices
            self.clustering.move_item(idx=oidx,cid=cid)
            self.G[REPRNAMES.CLUSTER_INTERACTION_INDICATOR].update_cnt_indicators()
            self.G.data.design_matrices[REPRNAMES.CLUSTER_INTERACTION_VALUE].M.build()
            if self._there_are_perturbations:
                self.G.data.design_matrices[REPRNAMES.PERT_VALUE].M.build()

            LOG_P.append(np.log(self.clustering.clusters[cid].size - 1) + \
                self.calculate_marginal_loglikelihood_slow_fast_sparse())
            LOG_KEYS.append(cid)


        # Calculate new cluster
        # =====================
        cid=self.clustering.make_new_cluster_with(idx=oidx)
        self.G[REPRNAMES.CLUSTER_INTERACTION_INDICATOR].update_cnt_indicators()
        self.G.data.design_matrices[REPRNAMES.CLUSTER_INTERACTION_VALUE].M.build()
        if self._there_are_perturbations:
            self.G.data.design_matrices[REPRNAMES.PERT_VALUE].M.build()
        
        LOG_KEYS.append(cid)
        LOG_P.append(np.log(concentration/self.m) + \
            self.calculate_marginal_loglikelihood_slow_fast_sparse())

        # Sample the assignment
        # =====================
        idx = sample_categorical_log(LOG_P)
        assigned_cid = LOG_KEYS[idx]
        curr_clus = self.clustering.idx2cid[oidx]

        if assigned_cid != curr_clus:
            self.clustering.move_item(idx=oidx,cid=assigned_cid)
            self.G[REPRNAMES.CLUSTER_INTERACTION_INDICATOR].update_cnt_indicators()

            # Change the mixing matrix for the interactions and (potentially) perturbations
            self.G.data.design_matrices[REPRNAMES.CLUSTER_INTERACTION_VALUE].M.build()
            if self._there_are_perturbations:
                self.G.data.design_matrices[REPRNAMES.PERT_VALUE].M.build()

    # @profile
    def calculate_marginal_loglikelihood_slow_fast(self):
        '''Marginalizes out the interactions and the perturbations
        '''
        # Build the parameters
        # ====================
        if self._there_are_perturbations:
            rhs = [REPRNAMES.PERT_VALUE, REPRNAMES.CLUSTER_INTERACTION_VALUE]
        else:
            rhs = [REPRNAMES.CLUSTER_INTERACTION_VALUE]
        
        
        y = self.y
        X = self.G.data.construct_rhs(keys=rhs, toarray=True)

        process_prec = self.process_prec
        prior_prec = build_prior_covariance(G=self.G, cov=False, order=rhs, sparse=False)
        prior_prec_diag = np.diag(prior_prec)
        prior_var = build_prior_covariance(G=self.G, cov=True, order=rhs, sparse=False)
        prior_mean = build_prior_mean(G=self.G, order=rhs, shape=(-1,1))

        # Calculate the marginalization
        # =============================
        a = X.T * process_prec

        beta_prec = a @ X + prior_prec
        beta_cov = pinv(beta_prec, self)
        beta_mean = beta_cov @ ( a @ y + prior_prec @ prior_mean )
        beta_mean = np.asarray(beta_mean).reshape(-1,1)

        try:
            beta_logdet = log_det(beta_cov, self)
        except:
            logging.critical('Crashed in log_det')
            logging.critical('beta_cov:\n{}'.format(beta_cov))
            logging.critical('prior_prec\n{}'.format(prior_prec))
            raise
        priorvar_logdet = log_det(prior_var, self)
        ll2 = 0.5 * (beta_logdet - priorvar_logdet)

        a = np.sum((prior_mean.ravel() ** 2) *prior_prec_diag)
        # np.asarray(prior_mean.T @ prior_prec @ prior_mean)[0,0]
        b = np.asarray(beta_mean.T @ beta_prec @ beta_mean)[0,0]
        ll3 = -0.5 * (a  - b)

        return ll2+ll3

    # @profile
    def calculate_marginal_loglikelihood_slow_fast_sparse(self):
        '''Marginalizes out the interactions and the perturbations
        '''
        # Build the parameters
        # ====================
        if self._there_are_perturbations:
            rhs = [REPRNAMES.PERT_VALUE, REPRNAMES.CLUSTER_INTERACTION_VALUE]
        else:
            rhs = [REPRNAMES.CLUSTER_INTERACTION_VALUE]
        
        
        y = self.y
        X = self.G.data.construct_rhs(keys=rhs, toarray=False)

        process_prec = self.process_prec_matrix
        prior_prec = build_prior_covariance(G=self.G, cov=False, order=rhs, sparse=True)  
        prior_prec_diag = build_prior_covariance(G=self.G, cov=False, order=rhs, diag=True)        
        prior_var = build_prior_covariance(G=self.G, cov=True, order=rhs, sparse=True)
        prior_mean = build_prior_mean(G=self.G, order=rhs, shape=(-1,1))

        # Calculate the marginalization
        # =============================

        # print('X')
        # print(type(X))
        # print(X.shape)

        # print('process_prec')
        # print(type(process_prec))
        # print(process_prec.shape)

        # print('prior_prec')
        # print(type(prior_prec))
        # print(prior_prec.shape)

        # print('prior mean')
        # print(type(prior_mean))
        # print(prior_mean.shape)

        a = X.T.dot(process_prec)
        beta_prec = a.dot(X) + prior_prec
        beta_cov = pinv(beta_prec, self)
        beta_mean = beta_cov @ ( a.dot(y) + prior_prec.dot(prior_mean))
        beta_mean = np.asarray(beta_mean).reshape(-1,1)

        try:
            beta_logdet = log_det(beta_cov, self)
        except:
            logging.critical('Crashed in log_det')
            logging.critical('beta_cov:\n{}'.format(beta_cov))
            logging.critical('prior_prec\n{}'.format(prior_prec))
            raise
        priorvar_logdet = log_det(prior_var, self)
        ll2 = 0.5 * (beta_logdet - priorvar_logdet)

        a = np.sum((prior_mean.ravel() ** 2) *prior_prec_diag)
        # np.asarray(prior_mean.T @ prior_prec @ prior_mean)[0,0]

        # print('beta_prec.shape', beta_prec.shape)
        # print('beta_mean.shape', beta_mean.shape)

        b = np.asarray(beta_mean.T @ beta_prec.dot(beta_mean))[0,0]
        ll3 = -0.5 * (a  - b)

        return ll2+ll3

    # Update MP - meant to be used during inference
    # =============================================
    def update_mp(self):
        '''Implements `update_slow` but parallelizes calculating the likelihood
        of being in a cluster. NOTE that this does not parallelize on the ASV level.

        On the first gibb step with initialize the workers that we implement with DASW (
        different arguments, single worker). For more information what this means look
        at pylab.multiprocessing documentation.

        If we initialized our pool as a pylab.multiprocessing.PersistentPool, then we 
        multiprocess the likelihood calculations for each asv. If we didnt then this 
        implementation has the same performance as `ClusterAssignments.update_slow_fast`.
        '''
        if self.G.data.zero_inflation_transition_policy is not None:
            raise NotImplementedError('Multiprocessing for zero inflation data is not implemented yet.' \
                ' Use `mp=None`')
        DMI = self.G.data.design_matrices[REPRNAMES.CLUSTER_INTERACTION_VALUE]
        DMP = self.G.data.design_matrices[REPRNAMES.PERT_VALUE]
        if self.clustering.n_clusters.sample_iter == 0 or self.pool == []:
            kwargs = {
                'n_asvs': len(self.G.data.asvs),
                'total_n_dts_per_asv': self.G.data.total_n_dts_per_asv,
                'n_replicates': self.G.data.n_replicates,
                'n_dts_for_replicate': self.G.data.n_dts_for_replicate,
                'there_are_perturbations': self._there_are_perturbations,
                'keypair2col_interactions': DMI.M.keypair2col,
                'keypair2col_perturbations': DMP.M.keypair2col,
                'n_perturbations': len(self.G.perturbations) if self._there_are_perturbations else None,
                'base_Xrows': DMI.base.rows,
                'base_Xcols': DMI.base.cols,
                'base_Xshape': DMI.base.shape,
                'base_Xpertrows': DMP.base.rows,
                'base_Xpertcols': DMP.base.cols,
                'base_Xpertshape': DMP.base.shape,
                'n_rowsM': DMI.M.n_rows,
                'n_rowsMpert': DMP.M.n_rows}

            if pl.ispersistentpool(self.pool):
                for _ in range(self.n_cpus):
                    self.pool.add_worker(SingleClusterFullParallelization(**kwargs))
            else:
                self.pool = SingleClusterFullParallelization(**kwargs)

        if self.clustering.n_clusters.sample_iter < self.delay:
            return

        if self.clustering.n_clusters.sample_iter % self.run_every_n_iterations != 0:
            return

        # Send in arguments for the start of the gibbs step
        start_time = time.time()
        base_Xdata = DMI.base.data
        self.concentration = self.G[REPRNAMES.CONCENTRATION].value
        y = self.G.data.construct_lhs(keys=[REPRNAMES.GROWTH_VALUE, REPRNAMES.SELF_INTERACTION_VALUE],
            kwargs_dict={REPRNAMES.GROWTH_VALUE: {'with_perturbations':False}})
        prior_var_interactions = self.G[REPRNAMES.PRIOR_VAR_INTERACTIONS].value
        prior_mean_interactions = self.G[REPRNAMES.PRIOR_MEAN_INTERACTIONS].value
        process_prec_diag = self.G[REPRNAMES.PROCESSVAR].prec

        if self._there_are_perturbations:
            prior_var_pert = self.G[REPRNAMES.PRIOR_VAR_PERT].get_single_value_of_perts()
            prior_mean_pert = self.G[REPRNAMES.PRIOR_MEAN_PERT].get_single_value_of_perts()
            base_Xpertdata = DMP.base.data
        else:
            prior_var_pert = None
            prior_mean_pert = None
            base_Xpertdata = None

        kwargs = {
            'base_Xdata': base_Xdata,
            'base_Xpertdata': base_Xpertdata,
            'concentration': self.concentration,
            'm': self.m,
            'y': y,
            'process_prec_diag': process_prec_diag,
            'prior_var_interactions': prior_var_interactions,
            'prior_var_pert': prior_var_pert,
            'prior_mean_interactions': prior_mean_interactions,
            'prior_mean_pert': prior_mean_pert}
        
        if pl.ispersistentpool(self.pool):
            self.pool.map('initialize_gibbs', [kwargs]*self.pool.num_workers)
        else:
            self.pool.initialize_gibbs(**kwargs)

        oidxs = npr.permutation(np.arange(len(self.G.data.asvs)))
        for iii, oidx in enumerate(oidxs):
            logging.info('{}/{} - {}'.format(iii, len(self.G.data.asvs), oidx))
            self.oidx = oidx
            self.gibbs_update_single_asv_parallel()

        self._strtime = time.time() - start_time

    def gibbs_update_single_asv_parallel(self):
        '''Update for a single asvs
        '''
        self.original_cluster = self.clustering.idx2cid[self.oidx]
        self.curr_cluster = self.original_cluster

        interactions = self.G[REPRNAMES.INTERACTIONS_OBJ]
        interaction_on_idxs = interactions.get_indicators(return_idxs=True)
        if self._there_are_perturbations:
            perturbation_on_idxs = [p.indicator.cluster_arg_array() for p in self.G.perturbations]
        else:
            perturbation_on_idxs = None

        # # Send topology parameters if the topology wont change
        # TODO: this works on windows but not on the cluster when dispatching?
        # if self.clustering.clusters[self.original_cluster].size > 1:
        #     use_saved_params = True
        #     kwargs = {
        #         'interaction_on_idxs': interaction_on_idxs,
        #         'perturbation_on_idxs': perturbation_on_idxs}
        #     if pl.ispersistentpool(self.pool):
        #         self.pool.map('initialize_oidx', [kwargs]*self.pool.num_workers)
        #     else:
        #         use_saved_params = False
        # else:
        #     use_saved_params = False
        use_saved_params = False

        if pl.ispersistentpool(self.pool):
            self.pool.staged_map_start('run')
        else:
            notpool_ret = []

        # Get the likelihood of the current configuration
        if self.clustering.clusters[self.original_cluster].size == 1:
            log_mult_factor = math.log(self.concentration/self.m)
        else:
            log_mult_factor = math.log(self.clustering.clusters[self.original_cluster].size - 1)

        if use_saved_params and pl.ispersistentpool(self.pool):
            interaction_on_idxs = None
            perturbation_on_idxs = None

        cluster_config = np.asarray([self.clustering.cid2cidx[self.clustering.idx2cid[i]] \
            for i in range(len(self.G.data.asvs))])

        kwargs = {
            'interaction_on_idxs': interaction_on_idxs,
            'perturbation_on_idxs': perturbation_on_idxs,
            'cluster_config': cluster_config,
            'log_mult_factor': log_mult_factor,
            'cid': self.original_cluster,
            'use_saved_params': use_saved_params}

        if pl.ispersistentpool(self.pool):
            self.pool.staged_map_put(kwargs)
        else:
            notpool_ret.append(self.pool.run(**kwargs))

        # Check every cluster
        for cid in self.clustering.order:
            if cid == self.original_cluster:
                continue
            self.clustering.move_item(idx=self.oidx, cid=cid)
            self.curr_cluster = cid

            if not use_saved_params:
                interaction_on_idxs = interactions.get_indicators(return_idxs=True)
                if self._there_are_perturbations:
                    perturbation_on_idxs = [p.indicator.cluster_arg_array() for p in self.G.perturbations]
                else:
                    perturbation_on_idxs = None

            cluster_config = np.asarray([self.clustering.cid2cidx[self.clustering.idx2cid[i]] \
                for i in range(len(self.G.data.asvs))])
            log_mult_factor = np.log(self.clustering.clusters[self.curr_cluster].size - 1)

            kwargs = {
                'interaction_on_idxs': interaction_on_idxs,
                'perturbation_on_idxs': perturbation_on_idxs,
                'cluster_config': cluster_config,
                'log_mult_factor': log_mult_factor,
                'cid': self.curr_cluster,
                'use_saved_params': use_saved_params}

            if pl.ispersistentpool(self.pool):
                self.pool.staged_map_put(kwargs)
            else:
                notpool_ret.append(self.pool.run(**kwargs))

        # Make a new cluster
        self.curr_cluster = self.clustering.make_new_cluster_with(idx=self.oidx)
        cluster_config = np.asarray([self.clustering.cid2cidx[self.clustering.idx2cid[i]] \
            for i in range(len(self.G.data.asvs))])
        interaction_on_idxs = interactions.get_indicators(return_idxs=True)
        if self._there_are_perturbations:
            perturbation_on_idxs = [p.indicator.cluster_arg_array() for p in self.G.perturbations]
        else:
            perturbation_on_idxs = None
        log_mult_factor = np.log(self.concentration/self.m)

        kwargs = {
            'interaction_on_idxs': interaction_on_idxs,
            'perturbation_on_idxs': perturbation_on_idxs,
            'cluster_config': cluster_config,
            'log_mult_factor': log_mult_factor,
            'cid': self.curr_cluster,
            'use_saved_params': False}

        # Put the values and get if necessary
        KEYS = []
        LOG_P = []
        if pl.ispersistentpool(self.pool):
            self.pool.staged_map_put(kwargs)
            ret = self.pool.staged_map_get()
        else:
            notpool_ret.append(self.pool.run(**kwargs))
            ret = notpool_ret
        for c, p in ret:
            KEYS.append(c)
            LOG_P.append(p)

        idx = sample_categorical_log(LOG_P)
        assigned_cid = KEYS[idx]

        if assigned_cid != self.original_cluster:
            logging.info('cluster changed')

        self.clustering.move_item(idx=self.oidx, cid=assigned_cid)

        self.G[REPRNAMES.CLUSTER_INTERACTION_INDICATOR].update_cnt_indicators()
        self.G.data.design_matrices[REPRNAMES.CLUSTER_INTERACTION_VALUE].M.build()
        if self._there_are_perturbations:
            self.G.data.design_matrices[REPRNAMES.PERT_VALUE].M.build()

        
class SingleClusterFullParallelization(pl.multiprocessing.PersistentWorker):
    '''Make the full parallelization
        - Mixture matricies for interactions and perturbations
        - calculating the marginalization for the sent in cluster

    Parameters
    ----------
    n_asvs : int
        Total number of OTUs
    total_n_dts_per_asv : int
        Total number of time changes for each OTU
    n_replicates : int
        Total number of replicates
    n_dts_for_replicate : np.ndarray
        Total number of time changes for each replicate
    there_are_perturbations : bool
        If True, there are perturbations
    keypair2col_interactions : np.ndarray
        These map the OTU indices of the pairs of OTUs to their column index
        in `big_X`
    keypair2col_perturbations : np.ndarray, None
        These map the OTU indices nad perturbation index to the column in
        `big_Xpert`. If there are no perturbations then this is None
    n_perturbations : int, None
        Number of perturbations. None if there are no perturbations
    base_Xrows, base_Xcols, base_Xpertrows, base_Xpertcols : np.ndarray
        These are the rows and columns necessary to build the interaction and perturbation 
        matrices, respectively. Whats passed in is the data vector and then we build
        it using sparse matrices
    n_rowsM, n_rowsMpert : int
        These are the number of rows for the mixing matrix for the interactions and
        perturbations respectively.
    '''
    def __init__(self, n_asvs, total_n_dts_per_asv, n_replicates, n_dts_for_replicate,
        there_are_perturbations, keypair2col_interactions, keypair2col_perturbations,
        n_perturbations, base_Xrows, base_Xcols, base_Xshape, base_Xpertrows, base_Xpertcols,
        base_Xpertshape, n_rowsM, n_rowsMpert):
        self.n_asvs = n_asvs
        self.total_n_dts_per_asv = total_n_dts_per_asv
        self.n_replicates = n_replicates
        self.n_dts_for_replicate = n_dts_for_replicate
        self.there_are_perturbations = there_are_perturbations
        self.keypair2col_interactions = keypair2col_interactions
        if self.there_are_perturbations:
            self.keypair2col_perturbations = keypair2col_perturbations
            self.n_perturbations = n_perturbations

        self.base_Xrows = base_Xrows
        self.base_Xcols = base_Xcols
        self.base_Xshape = base_Xshape
        self.base_Xpertrows = base_Xpertrows
        self.base_Xpertcols = base_Xpertcols
        self.base_Xpertshape = base_Xpertshape

        self.n_rowsM = n_rowsM
        self.n_rowsMpert = n_rowsMpert

    def initialize_gibbs(self, base_Xdata, base_Xpertdata, concentration, m, y,
        process_prec_diag, prior_var_interactions, prior_var_pert, 
        prior_mean_interactions, prior_mean_pert):
        '''Pass in the information that changes every Gibbs step

        Parameters
        ----------
        base_X : scipy.sparse.csc_matrix
            Sparse matrix for the interaction terms
        base_Xpert : scipy.sparse.csc_matrix, None
            Sparse matrix for the perturbation terms
            If None, there are no perturbations
        concentration : float
            This is the concentration of the system
        m : int
            This is the auxiliary variable for the marginalization
        y : np.ndarray
            This is the observation array
        process_prec_diag : np.ndarray
            This is the process precision diagonal
        '''
        self.base_X = scipy.sparse.coo_matrix(
            (base_Xdata,(self.base_Xrows,self.base_Xcols)),
            shape=self.base_Xshape).tocsc()
        
        self.concentration = concentration
        self.m = m
        self.y = y.reshape(-1,1)
        self.n_rows = len(y)
        self.n_cols_X = self.base_X.shape[1]
        self.prior_var_interactions = prior_var_interactions
        self.prior_prec_interactions = 1/prior_var_interactions
        self.prior_mean_interactions = prior_mean_interactions

        if self.there_are_perturbations:
            self.base_Xpert = scipy.sparse.coo_matrix(
                (base_Xpertdata,(self.base_Xpertrows,self.base_Xpertcols)),
                shape=self.base_Xpertshape).tocsc()
            self.prior_var_pert = prior_var_pert
            self.prior_prec_pert = 1/prior_var_pert
            self.prior_mean_pert = prior_mean_pert

        self.process_prec_matrix = scipy.sparse.dia_matrix(
            (process_prec_diag,[0]), shape=(len(process_prec_diag),len(process_prec_diag))).tocsc()

    def initialize_oidx(self, interaction_on_idxs, perturbation_on_idxs):
        '''Pass in the parameters that change for every OTU - potentially

        Parameters
        ----------
        
        '''
        self.saved_interaction_on_idxs = interaction_on_idxs
        self.saved_perturbation_on_idxs = perturbation_on_idxs

    # @profile
    def run(self, interaction_on_idxs, perturbation_on_idxs, cluster_config, log_mult_factor, cid, 
        use_saved_params):
        '''Pass in the parameters for the specific cluster assignment for
        the OTU and run the marginalization


        Parameters
        ----------
        interaction_on_idxs : np.array(int)
            An array of indices for the interactions that are on. Assumes that the
            clustering is in the order specified in `cluster_config`
        perturbation_on_idxs : list(np.ndarray(int)), None
            If there are perturbations, then we set the perturbation idxs on
            Each element in the list are the indices of that perturbation that are on
        cluster_config : list(list(int))
            This is the cluster configuration and in cluster order.
        log_mult_factor : float
            This is the log multiplication factor that we add onto the marginalization
        use_saved_params : bool
            If True, passed in `interaction_on_idxs` and `perturbation_on_idxs` are None
            and we can use `saved_interaction_on_idxs` and `saved_perturbation_on_idxs`
        '''
        if use_saved_params:
            interaction_on_idxs = self.saved_interaction_on_idxs
            perturbation_on_idxs = self.saved_perturbation_on_idxs

        # We need to make the arrays for interactions and perturbations
        self.set_clustering(cluster_config=cluster_config)
        Xinteractions = self.build_interactions_matrix(on_columns=interaction_on_idxs)

        if self.there_are_perturbations:
            Xperturbations = self.build_perturbations_matrix(on_columns=perturbation_on_idxs)
            X = scipy.sparse.hstack([Xperturbations, Xinteractions])
        else:
            X = Xinteractions
        self.X = X
        self.prior_mean = self.build_prior_mean(on_interactions=interaction_on_idxs,
            on_perturbations=perturbation_on_idxs)
        self.prior_cov, self.prior_prec, self.prior_prec_diag = self.build_prior_cov_and_prec_and_diag(
            on_interactions=interaction_on_idxs, on_perturbations=perturbation_on_idxs)
        
        return cid, self.calculate_marginal_loglikelihood_slow_fast_sparse()

    def set_clustering(self, cluster_config):
        self.clustering = CondensedClustering(oidx2cidx=cluster_config)
        self.iidx2cidxpair = np.zeros(shape=(len(self.clustering)*(len(self.clustering)-1), 2), 
            dtype=int)
        self.iidx2cidxpair = SingleClusterFullParallelization.make_iidx2cidxpair(
            ret=self.iidx2cidxpair,
            n_clusters=len(self.clustering))

    def build_interactions_matrix(self, on_columns):
        '''Build the interaction matrix

        First we make the rows and columns for the mixing matrix,
        then we multiple the base matrix and the mixing matrix.
        '''
        rows = []
        cols = []

        # c2ciidx = Cluster-to-Cluster Interaction InDeX
        c2ciidx = 0
        for ccc in on_columns:
            tcidx = self.iidx2cidxpair[ccc, 0]
            scidx = self.iidx2cidxpair[ccc, 1]
            
            smems = self.clustering.clusters[scidx]
            tmems = self.clustering.clusters[tcidx]
            
            a = np.zeros(len(smems)*len(tmems), dtype=int)
            rows.append(SingleClusterFullParallelization.get_indices(a,
                self.keypair2col_interactions, tmems, smems))
            cols.append(np.full(len(tmems)*len(smems), fill_value=c2ciidx))
            c2ciidx += 1

        rows = np.asarray(list(itertools.chain.from_iterable(rows)))
        cols = np.asarray(list(itertools.chain.from_iterable(cols)))
        data = np.ones(len(rows), dtype=int)

        # print('rows', rows)
        # print(cols)
        # print(data)
        # print((self.n_rowsM, c2ciidx))

        M = scipy.sparse.coo_matrix((data,(rows,cols)),
            shape=(self.n_rowsM, c2ciidx)).tocsc()
        ret = self.base_X @ M
        return ret

    def build_perturbations_matrix(self, on_columns):
        if not self.there_are_perturbations:
            raise ValueError('You should not be here')
        
        keypair2col = self.keypair2col_perturbations
        rows = []
        cols = []

        col = 0
        for pidx, pert_ind_idxs in enumerate(on_columns):
            for cidx in pert_ind_idxs:
                for oidx in self.clustering.clusters[cidx]:
                    rows.append(keypair2col[oidx, pidx])
                    cols.append(col)
                col += 1
        
        data = np.ones(len(rows), dtype=np.float64)
        M = scipy.sparse.coo_matrix((data,(rows,cols)),
            shape=(self.n_rowsMpert, col)).tocsc()
        ret = self.base_Xpert @ M
        return ret

    def build_prior_mean(self, on_interactions, on_perturbations):
        '''Build the prior mean array

        Perturbations go first and then interactions
        '''
        ret = []
        for pidx, pert in enumerate(on_perturbations):
            ret = np.append(ret, 
                np.full(len(pert), fill_value=self.prior_mean_pert[pidx]))
        ret = np.append(ret, np.full(len(on_interactions), fill_value=self.prior_mean_interactions))
        return ret.reshape(-1,1)

    def build_prior_cov_and_prec_and_diag(self, on_interactions, on_perturbations):
        '''Build the prior covariance matrices and others
        '''
        ret = []
        for pidx, pert in enumerate(on_perturbations):
            ret = np.append(ret, 
                np.full(len(pert), fill_value=self.prior_var_pert[pidx]))
        prior_var_diag = np.append(ret, np.full(len(on_interactions), fill_value=self.prior_var_interactions))
        prior_prec_diag = 1/prior_var_diag

        prior_var = scipy.sparse.dia_matrix((prior_var_diag,[0]), 
            shape=(len(prior_var_diag),len(prior_var_diag))).tocsc()
        prior_prec = scipy.sparse.dia_matrix((prior_prec_diag,[0]), 
            shape=(len(prior_prec_diag),len(prior_prec_diag))).tocsc()

        return prior_var, prior_prec, prior_prec_diag

    # @profile
    def calculate_marginal_loglikelihood_slow_fast_sparse(self):
        y = self.y
        X = self.X
        process_prec = self.process_prec_matrix
        prior_mean = self.prior_mean
        prior_cov = self.prior_cov
        prior_prec = self.prior_prec
        prior_prec_diag = self.prior_prec_diag

        a = X.T.dot(process_prec)
        beta_prec = a.dot(X) + prior_prec
        beta_cov = pinv(beta_prec, self)
        beta_mean = beta_cov @ (a.dot(y) + prior_prec.dot(prior_mean))
        beta_mean = np.asarray(beta_mean).reshape(-1,1)

        try:
            beta_logdet = log_det(beta_cov, self)
        except:
            logging.critical('Crashed in log_det')
            logging.critical('beta_cov:\n{}'.format(beta_cov))
            logging.critical('prior_prec\n{}'.format(prior_prec))
            raise
        
        priorvar_logdet = log_det(prior_cov, self)
        ll2 = 0.5 * (beta_logdet - priorvar_logdet)

        bEbprior = np.asarray(prior_mean.T @ prior_prec.dot(prior_mean))[0,0]
        bEb = np.asarray(beta_mean.T @ beta_prec.dot(beta_mean) )[0,0]
        ll3 = 0.5 * (bEb - bEbprior)

        self.a = a
        self.beta_prec = beta_prec

        return ll2 + ll3

    @staticmethod
    @numba.jit(nopython=True, cache=True)
    def get_indices(a, keypair2col, tmems, smems):
        '''Use Just in Time compilation to reduce the 'getting' time
        by about 95%

        Parameters
        ----------
        keypair2col : np.ndarray
            Maps (target_oidx, source_oidx) pair to the place the interaction
            index would be on a full interactio design matrix on the OTU level
        tmems, smems : np.ndarray
            These are the OTU indices in the target cluster and the source cluster
            respectively
        '''
        i = 0
        for tidx in tmems:
            for sidx in smems:
                a[i] = keypair2col[tidx, sidx]
                i += 1
        return a

    @staticmethod
    # @numba.jit(nopython=True, cache=True)
    def make_iidx2cidxpair(ret, n_clusters):
        '''Map the index of a cluster interaction to (dst,src) of clusters

        Parameters
        ----------
        n_clusters : int
            Number of clusters

        Returns
        -------
        np.ndarray(n_interactions,2)
            First column is the destination cluster index, second column is the source
            cluster index
        '''
        i = 0
        for dst_cidx in range(n_clusters):
            for src_cidx in range(n_clusters):
                if dst_cidx == src_cidx:
                    continue
                ret[i,0] = dst_cidx
                ret[i,1] = src_cidx
                i += 1
        
        return ret


class CondensedClustering:
    '''Condensed clustering object that is not associated with the graph

    Parameters
    ----------
    oidx2cidx : np.ndarray
        Maps the cluster assignment to each asv.
        index -> ASV index
        output -> cluster index

    '''
    def __init__(self, oidx2cidx):

        self.clusters = []
        self.oidx2cidx = oidx2cidx
        a = {}
        for oidx, cidx in enumerate(self.oidx2cidx):
            if cidx not in a:
                a[cidx] = [oidx]
            else:
                a[cidx].append(oidx)
        
        cidx = 0
        while cidx in a:
            self.clusters.append(np.asarray(a[cidx], dtype=int))
            cidx += 1

    def __len__(self):
        return len(self.clusters)


# Filtering
# ---------
class TrajectorySet(pl.graph.Node):
    '''This aggregates a set of trajectories from each set

    Parameters
    ----------
    name : str
        Name of the object
    G : pylab.graph.Graph
        Graph object to attach it to
    '''
    def __init__(self, name, G, **kwargs):
        pl.graph.Node.__init__(self, name=name, G=G)
        self.value = []
        n_asvs = self.G.data.n_asvs

        for ridx, subj in enumerate(self.G.data.subjects):
            n_timepoints = self.G.data.n_timepoints_for_replicate[ridx]

            # initialize values to zeros for initialization
            self.value.append(pl.variables.Variable(
                name=name+'_{}'.format(subj.name), G=G, shape=(n_asvs, n_timepoints),
                value=np.zeros((n_asvs, n_timepoints), dtype=float), **kwargs))
        prior = pl.variables.Normal(
            mean=pl.variables.Constant(name=self.name+'_prior_mean', value=0, G=self.G),
            var=pl.variables.Constant(name=self.name+'_prior_var', value=1, G=self.G),
            name=self.name+'_prior', G=self.G)
        self.add_prior(prior)

    def __getitem__(self, ridx):
        return self.value[ridx]

    @property
    def sample_iter(self):
        return self.value[0].sample_iter

    def reset_value_size(self):
        '''Change the size of the trajectory when we set the intermediate timepoints
        '''
        n_asvs = self.G.data.n_asvs
        for ridx in range(len(self.value)):
            n_timepoints = self.G.data.n_timepoints_for_replicate[ridx]
            self.value[ridx].value = np.zeros((n_asvs, n_timepoints),dtype=float)
            self.value[ridx].set_value_shape(self.value[ridx].value.shape)

    def _vectorize(self):
        '''Get all the data in vector form
        '''
        vals = np.array([])
        for data in self.value:
            vals = np.append(vals, data.value)
        return vals

    def mean(self):
        return np.mean(self._vectorize())

    def var(self):
        return np.var(self._vectorize())

    def iter_indices(self):
        '''Iterate through the indices and the values
        '''
        for ridx in range(self.G.data.n_replicates):
            for tidx in range(self.G.data.n_timepoints_for_replicate[ridx]):
                for oidx in range(self.G.data.asvs.n_asvs):
                    yield (ridx, tidx, oidx)

    def set_trace(self, *args, **kwargs):
        for ridx in range(len(self.value)):
            self.value[ridx].set_trace(*args, **kwargs)

    def add_trace(self):
        for ridx in range(len(self.value)):
            # Set the zero inflation values to nans
            self.value[ridx].value[~self.G[REPRNAMES.ZERO_INFLATION].value[ridx]] = np.nan
            self.value[ridx].add_trace()


class FilteringLogMP(pl.graph.Node):
    '''This is the posterior for the latent trajectory that are
    sampled using a standard normal Metropolis-Hastings proposal.

    This is the multiprocessing version of the class. All of the computation is
    done on the subject level in parallel.

    Parallelization Modes
    ---------------------
    'debug'
        If this is selected, then we dont actually parallelize, but we go in
        order of the objects in sequential order. We would do this if we want
        to benchmark within each processor or do easier print statements
    'full'
        This is where each subject gets their own process

    This assumes that we are using the log model
    '''
    def __init__(self, mp, zero_inflation_transition_policy,**kwargs):
        '''
        Parameters
        ----------
        mp : str
            'debug'
                Does not actually parallelize, does it serially - we do this in case we
                want to debug and/or benchmark
            'full'
                Send each replicate to a processor each
        zero_inflation_transition_policy : None, str
            Type of zero inflation to do. If None then there is no zero inflation
        '''
        kwargs['name'] = STRNAMES.FILTERING
        pl.graph.Node.__init__(self, **kwargs)
        self.x = TrajectorySet(name=STRNAMES.LATENT_TRAJECTORY, G=self.G)
        self.mp = mp
        self.zero_inflation_transition_policy = zero_inflation_transition_policy

        self.print_vals = False
        self._strr = 'parallel'

    def __str__(self):
        return self._strr

    @property
    def sample_iter(self):
        # It doesnt matter if we chose q or x because they are both the same
        return self.x.sample_iter

    def initialize(self, x_value_option, a0, a1, v1, v2, essential_timepoints, tune, 
        proposal_init_scale, intermediate_step, intermediate_interpolation=None, 
        delay=0, bandwidth=None, window=None, target_acceptance_rate=0.44, 
        calculate_qpcr_loglik=True):
        '''Initialize the values of the error model (values for the
        latent and the auxiliary trajectory). Additionally this sets
        the intermediate time points

        Initialize the values of the prior.

        Parameters
        ----------
        x_value_option : str
            Option to initialize the value of the latent trajectory.
            Options
                'coupling'
                    Sample the values around the data with extremely low variance.
                    This also truncates the data so that it stays > 0.
                'moving-avg'
                    Initialize the values using a moving average around the points.
                    The bandwidth of the filter is by number of days, not the order
                    of timepoints. You must also provide the argument `bandwidth`.
                'loess', 'auto'
                    Implements the initialization of the values using LOESS (Locally
                    Estimated Scatterplot Smoothing) algorithm. You must also provide
                    the `window` parameter
        tune : tuple(int, int)
            This is how often to tune the individual covariances
            The first element indicates which MCMC sample to stop the tuning
            The second element is how often to update the proposal covariance
        a0, a1 : float, str
            These are the hyperparameters to calculate the dispersion of the
            negative binomial.
        v1, v2 : float, int, str
            These are the values used to calulcate the coupling variance between
            x and q
        intermediate_step : tuple(str, args), array, None
            This is the type of interemediate timestep to intialize and the arguments
            for them. If this is None, then we do no intermediate timesteps.
            Options:
                'step'
                    args: (stride (numeric), eps (numeric))
                    We simulate at each timepoint every `stride` days.
                    We do not set an intermediate time point if it is within `eps`
                    days of a given data point.
                'preserve-density'
                    args: (n (int), eps (numeric))
                    We preserve the denisty of the given data by only simulating data
                    `n` times between each essential datapoint. If a timepoint is within
                    `eps` days of a given timepoint then we do not make an intermediate
                    point there.
                'manual;
                    args: np.ndarray
                    These are the points that we want to set. If these are not given times
                    then we set them as timepoints
        intermediate_interpolation : str
            This is the type of interpolation to perform on the intermediate timepoints.
            Options:
                'linear-interpolation', 'auto'
                    Perform linear interpolation between the two closest given timepoints
        essential_timepoints : np.ndarray, str, None
            These are the timepoints that must be included in each subject. If one of the
            subjects has a missing timepoint there then we use an intermediate time point
            that this timepoint. It is initialized with linear interpolation. If all of the
            timepoints specified in this vector are included in a subject then nothing is
            done. If it is a str:
                'union', 'auto'
                    We take a union of all the timepoints in each subject and make sure
                    that all of the subjects have all those points.
        bandwidth : float
            This is the day bandwidth of the filter if the initialization method is
            done with 'moving-avg'.
        window : int
            This is the window term for the LOESS initialization scheme. This is
            only used if value_initialization is done with 'loess'
        target_acceptance_rate : numeric
            This is the target acceptance rate for each time point individually
        calculate_qpcr_loglik : bool
            If True, calculate the loglikelihood of the qPCR measurements during the 
            proposal
        '''
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        self.delay = delay
        self._there_are_perturbations = self.G.perturbations is not None

        # Set the hyperparameters
        if not pl.isfloat(target_acceptance_rate):
            raise TypeError('`target_acceptance_rate` must be a float'.format(
                type(target_acceptance_rate)))
        if target_acceptance_rate < 0 or target_acceptance_rate > 1:
            raise ValueError('`target_acceptance_rate` ({}) must be in (0,1)'.format(
                target_acceptance_rate))
        if not pl.istuple(tune):
            raise TypeError('`tune` ({}) must be a tuple'.format(type(tune)))
        if len(tune) != 2:
            raise ValueError('`tune` ({}) must have 2 elements'.format(len(tune)))
        if not pl.isint(tune[0]):
            raise TypeError('`tune` ({}) 1st parameter must be an int'.format(type(tune[0])))
        if tune[0] < 0:
            raise ValueError('`tune` ({}) 1st parameter must be > 0'.format(tune[0]))
        if not pl.isint(tune[1]):
            raise TypeError('`tune` ({}) 2nd parameter must be an int'.format(type(tune[1])))
        if tune[1] < 0:
            raise ValueError('`tune` ({}) 2nd parameter must be > 0'.format(tune[1]))
        
        if not pl.isnumeric(a0):
            raise TypeError('`a0` ({}) must be a numeric type'.format(type(a0)))
        elif a0 <= 0:
            raise ValueError('`a0` ({}) must be > 0'.format(a0))
        if not pl.isnumeric(a1):
            raise TypeError('`a1` ({}) must be a numeric type'.format(type(a1)))
        elif a1 <= 0:
            raise ValueError('`a1` ({}) must be > 0'.format(a1))

        if not pl.isnumeric(proposal_init_scale):
            raise TypeError('`proposal_init_scale` ({}) must be a numeric type (int, float)'.format(
                type(proposal_init_scale)))
        if proposal_init_scale < 0:
            raise ValueError('`proposal_init_scale` ({}) must be positive'.format(
                proposal_init_scale))

        self.tune = tune
        self.a0 = a0
        self.a1 = a1
        self.target_acceptance_rate = target_acceptance_rate
        self.proposal_init_scale = proposal_init_scale
        self.v1 = v1
        self.v2 = v2

        # Set the essential timepoints (check to see if there is any missing data)
        if essential_timepoints is not None:
            logging.info('Setting up the essential timepoints')
            if pl.isstr(essential_timepoints):
                if essential_timepoints in ['auto', 'union']:
                    essential_timepoints = set()
                    for ts in self.G.data.times:
                        essential_timepoints = essential_timepoints.union(set(list(ts)))
                    essential_timepoints = np.sort(list(essential_timepoints))
                else:
                    raise ValueError('`essential_timepoints` ({}) not recognized'.format(
                        essential_timepoints))
            elif not pl.isarray(essential_timepoints):
                raise TypeError('`essential_timepoints` ({}) must be a str or an array'.format(
                    type(essential_timepoints)))
            logging.info('Essential timepoints: {}'.format(essential_timepoints))
            self.G.data.set_timepoints(times=essential_timepoints, eps=None, reset_timepoints=True)
            self.x.reset_value_size()

        # Set the intermediate timepoints if necessary
        if intermediate_step is not None:
            # Set the intermediate timepoints in the data
            if not pl.istuple(intermediate_step):
                raise TypeError('`intermediate_step` ({}) must be a tuple'.format(
                    type(intermediate_step)))
            if len(intermediate_step) != 2:
                raise ValueError('`intermediate_step` ({}) must be length 2'.format(
                    len(intermediate_step)))
            f, args = intermediate_step
            if not pl.isstr(f):
                raise TypeError('intermediate_step type ({}) must be a str'.format(type(f)))
            if f == 'step':
                if not pl.istuple(args):
                    raise TypeError('`args` ({}) must be a tuple'.format(type(args)))
                if len(args) != 2:
                    raise TypeError('`args` ({}) must have 2 arguments'.format(len(args)))
                step, eps = args
                self.G.data.set_timepoints(timestep=step, eps=eps, reset_timepoints=False)
            elif f == 'preserve-density':
                if not pl.istuple(args):
                    raise TypeError('`args` ({}) must be a tuple'.format(type(args)))
                if len(args) != 2:
                    raise TypeError('`args` ({}) must have 2 arguments'.format(len(args)))
                n, eps = args
                if not pl.isint(n):
                    raise TypeError('`n` ({}) must be an int'.format(type(n)))

                # For each timepoint, add `n` intermediate timepoints
                for ridx in range(self.G.data.n_replicates):
                    times = []
                    for i in range(len(self.G.data.times[ridx])-1):
                        t0 = self.G.data.times[ridx][i]
                        t1 = self.G.data.times[ridx][i+1]
                        step = (t1-t0)/(n+1)
                        times = np.append(times, np.arange(t0,t1,step=step))
                    times = np.sort(np.unique(times))
                    # print('\n\ntimes to put in', times)
                    self.G.data.set_timepoints(times=times, eps=eps, ridx=ridx, reset_timepoints=False)
                    # print('times for ridx {}'.format(self.G.data.times[ridx]))
                    # print('len times', len(self.G.data.times[ridx]))
                    # print('data shape', self.G.data.data[ridx].shape)

                # sys.exit()
            elif f == 'manual':
                raise NotImplementedError('Not Implemented')
            else:
                raise ValueError('`intermediate_step type ({}) not recognized'.format(f))
            self.x.reset_value_size()

        if intermediate_interpolation is not None:
            if intermediate_interpolation in ['linear-interpolation', 'auto']:
                for ridx in range(self.G.data.n_replicates):
                    for tidx in range(self.G.data.n_timepoints_for_replicate[ridx]):
                        if tidx not in self.G.data.given_timeindices[ridx]:
                            # We need to interpolate this time point
                            # get the previous given and next given timepoint
                            prev_tidx = None
                            for ii in range(tidx-1,-1,-1):
                                if ii in self.G.data.given_timeindices[ridx]:
                                    prev_tidx = ii
                                    break
                            if prev_tidx is None:
                                # Set to the same as the closest forward timepoint then continue
                                next_idx = None
                                for ii in range(tidx+1, self.G.data.n_timepoints_for_replicate[ridx]):
                                    if ii in self.G.data.given_timeindices[ridx]:
                                        next_idx = ii
                                        break
                                self.G.data.data[ridx][:,tidx] = self.G.data.data[ridx][:,next_idx]
                                continue

                            next_tidx = None
                            for ii in range(tidx+1, self.G.data.n_timepoints_for_replicate[ridx]):
                                if ii in self.G.data.given_timeindices[ridx]:
                                    next_tidx = ii
                                    break
                            if next_tidx is None:
                                # Set to the previous timepoint then continue
                                self.G.data.data[ridx][:,tidx] = self.G.data.data[ridx][:,prev_tidx]
                                continue

                            # Interpolate from prev_tidx to next_tidx
                            x = self.G.data.times[ridx][tidx]
                            x0 = self.G.data.times[ridx][prev_tidx]
                            y0 = self.G.data.data[ridx][:,prev_tidx]
                            x1 = self.G.data.times[ridx][next_tidx]
                            y1 = self.G.data.data[ridx][:,next_tidx]
                            self.G.data.data[ridx][:,tidx] = y0 * (1-((x-x0)/(x1-x0))) + y1 * (1-((x1-x)/(x1-x0)))
            else:
                raise ValueError('`intermediate_interpolation` ({}) not recognized'.format(intermediate_interpolation))

        # Initialize the latent trajectory
        if not pl.isstr(x_value_option):
            raise TypeError('`x_value_option` ({}) is not a str'.format(type(x_value_option)))
        if x_value_option == 'coupling':
            self._init_coupling()
        elif x_value_option == 'moving-avg':
            if not pl.isnumeric(bandwidth):
                raise TypeError('`bandwidth` ({}) must be a numeric'.format(type(bandwidth)))
            if bandwidth <= 0:
                raise ValueError('`bandwidth` ({}) must be positive'.format(bandwidth))
            self.bandwidth = bandwidth
            self._init_moving_avg()
        elif x_value_option in ['loess', 'auto']:
            if window is None:
                raise TypeError('If `value_option` is loess, then `window` must be specified')
            if not pl.isint(window):
                raise TypeError('`window` ({}) must be an int'.format(type(window)))
            if window <= 0:
                raise ValueError('`window` ({}) must be > 0'.format(window))
            self.window = window
            self._init_loess()
        else:
            raise ValueError('`x_value_option` ({}) not recognized'.format(x_value_option))

        # Get necessary data and set the parallel objects
        if self._there_are_perturbations:
            pert_starts = []
            pert_ends = []
            for perturbation in self.G.perturbations:
                pert_starts.append(perturbation.starts)
                pert_ends.append(perturbation.ends)
        else:
            pert_starts = None
            pert_ends = None

        if self.mp is None:
            self.mp = 'debug'
        if not pl.isstr(self.mp):
            raise TypeError('`mp` ({}) must either be a string or None'.format(type(self.mp)))
        if self.mp == 'debug':
            self.pool = []
        elif self.mp == 'full':
            self.pool = pl.multiprocessing.PersistentPool(G=self.G, ptype='sadw')
            self.worker_pids = []
        else:
            raise ValueError('`mp` ({}) not recognized'.format(self.mp))

        for ridx, subj in enumerate(self.G.data.subjects):
            # Set up qPCR measurements and reads to send
            qpcr_log_measurements = {}
            for t in self.G.data.given_timepoints[ridx]:
                qpcr_log_measurements[t] = self.G.data.qpcr[ridx][t].log_data
            reads = self.G.data.subjects.iloc(ridx).reads

            worker = SubjectLogTrajectorySetMP()
            worker.initialize(
                zero_inflation_transition_policy=self.zero_inflation_transition_policy,
                times=self.G.data.times[ridx],
                qpcr_log_measurements=qpcr_log_measurements,
                reads=reads,
                there_are_intermediate_timepoints=True,
                there_are_perturbations=self._there_are_perturbations,
                pv_global=self.G[REPRNAMES.PROCESSVAR].global_variance,
                x_prior_mean=np.log(1e7),
                x_prior_std=1e10,
                tune=tune[1],
                delay=delay,
                end_iter=tune[0],
                proposal_init_scale=proposal_init_scale,
                a0=a0,
                a1=a1,
                x=self.x[ridx].value,
                pert_starts=np.asarray(pert_starts),
                pert_ends=np.asarray(pert_ends),
                ridx=ridx,
                subjname=subj.name,
                calculate_qpcr_loglik=calculate_qpcr_loglik,
                h5py_xname=self.x[ridx].name,
                target_acceptance_rate=self.target_acceptance_rate)
            if self.mp == 'debug':
                self.pool.append(worker)
            elif self.mp == 'full':
                pid = self.pool.add_worker(worker)
                self.worker_pids.append(pid)

        # Set the data to the latent values
        self.set_latent_as_data(update_values=False)

        self.total_n_datapoints = 0
        for ridx in range(self.G.data.n_replicates):
            self.total_n_datapoints += self.x[ridx].value.shape[0] * self.x[ridx].value.shape[1]

    def _init_coupling(self):
        '''Initialize `x` by sampling around the data using a small
        variance using a truncated normal distribution
        '''
        for ridx, tidx, oidx in self.x.iter_indices():
            val = self.G.data.data[ridx][oidx,tidx]
            self.x[ridx][oidx,tidx] = pl.random.truncnormal.sample(
                mean=val,
                std=math.sqrt(self.v1 * (val ** 2) + self.v2),
                low=0, high=float('inf'))

    def _init_moving_avg(self):
        '''Initializes `x` by using a moving
        average over the data - using `self.bandwidth` as the bandwidth
        of number of days - it then samples around that point using the
        coupling variance.

        If there are no other points within the bandwidth around the point,
        then it just samples around the current timepoint with the coupling
        variance.
        '''
        for ridx in range(self.G.data.n_replicates):
            for tidx in range(self.G.data.n_timepoints_for_replicate[ridx]):
                tidx_low = np.searchsorted(
                    self.G.data.times[ridx], self.G.data.times[ridx][tidx]-self.bandwidth)
                tidx_high = np.searchsorted(
                    self.G.data.times[ridx], self.G.data.times[ridx][tidx]+self.bandwidth)

                for oidx in range(len(self.G.data.asvs)):
                    val = np.mean(self.G.data.data[ridx][oidx, tidx_low: tidx_high])
                    self.x[ridx][oidx,tidx] = pl.random.truncnormal.sample(
                        mean=val,
                        std=math.sqrt(self.v1 * (val ** 2) + self.v2),
                        low=0, high=float('inf'))

    def _init_loess(self):
        '''Initialize the data using LOESS algorithm and then samples around that
        the coupling variance we implement the LOESS algorithm in the module
        `fit_loess.py`
        '''
        for ridx in range(self.G.data.n_replicates):
            xx = self.G.data.times[ridx]
            for oidx in range(len(self.G.data.asvs)):
                yy = self.G.data.data[ridx][oidx, :]
                loess = _Loess(xx, yy)

                for tidx, t in enumerate(self.G.data.times[ridx]):
                    val = loess.estimate(t, window=self.window)
                    self.x[ridx][oidx,tidx] = pl.random.truncnormal.sample(
                        mean=val,
                        std=math.sqrt(self.v1 * (val ** 2) + self.v2),
                        low=0, high=float('inf'))

                    if np.isnan(self.x[ridx][oidx, tidx]):
                        print('crashed here', ridx, tidx, oidx)
                        print('mean', val)
                        print('t', t)
                        print('yy', yy)
                        print('std', math.sqrt(self.v1 * (val ** 2) + self.v2))
                        raise ValueError('')

    def set_latent_as_data(self, update_values=True):
        '''Change the values in the data matrix so that it is the latent variables
        '''
        data = []
        for obj in self.x.value:
            data.append(obj.value)
        self.G.data.data = data
        if update_values:
            self.G.data.update_values()

    def add_trace(self):
        self.x.add_trace()

    def set_trace(self, *args, **kwargs):
        self.x.set_trace(*args, **kwargs)

    def kill(self):
        if self.mp == 'full':
            self.pool.kill()

    def update(self):
        '''Send out to each parallel object
        '''
        if self.sample_iter < self.delay:
            return
        start_time = time.time()

        growth = self.G[REPRNAMES.GROWTH_VALUE].value.ravel()
        self_interactions = self.G[REPRNAMES.SELF_INTERACTION_VALUE].value.ravel()
        pv = self.G[REPRNAMES.PROCESSVAR].value
        interactions = self.G[REPRNAMES.INTERACTIONS_OBJ].get_datalevel_value_matrix(
            set_neg_indicators_to_nan=False)
        perts = None
        if self._there_are_perturbations:
            perts = []
            for perturbation in self.G.perturbations:
                perts.append(perturbation.item_array().reshape(-1,1))
            perts = np.hstack(perts)

        # zero_inflation = [self.G[REPRNAMES.ZERO_INFLATION].value[ridx] for ridx in range(self.G.data.n_replicates)]
        qpcr_vars = []
        for aaa in self.G[REPRNAMES.QPCR_VARIANCES].value:
            qpcr_vars.append(aaa.value)
        
        
        kwargs = {'growth':growth, 'self_interactions':self_interactions,
            'pv':pv, 'interactions':interactions, 'perturbations':perts, 
            'zero_inflation_data': None, 'qpcr_variances':qpcr_vars}

        str_acc = [None]*self.G.data.n_replicates
        if self.mp == 'debug':

            for ridx in range(self.G.data.n_replicates):
                _, x, acc_rate = self.pool[ridx].persistent_run(**kwargs)
                self.x[ridx].value = x
                str_acc[ridx] = '{:.3f}'.format(acc_rate)

        else:
            # raise NotImplementedError('Multiprocessing for filtering with zero inflation ' \
            #     'is not implemented')
            ret = self.pool.map(func='persistent_run', args=kwargs)
            for ridx, x, acc_rate in ret:
                self.x[ridx].value = x
                str_acc[ridx] = '{:.3f}'.format(acc_rate)

        self.set_latent_as_data()

        t = time.time() - start_time
        try:
            self._strr = 'Time: {:.4f}, Acc: {}, data/sec: {:.2f}'.format(t,
                str(str_acc).replace("'",''), self.total_n_datapoints/t)
        except:
            self._strr = 'NA'


class SubjectLogTrajectorySetMP(pl.multiprocessing.PersistentWorker):
    '''This performs filtering on a multiprocessing level. We send the
    other parameters of the model and return the filtered `x` and values.
    With multiprocessing, this class has ~91% efficiency. Additionally, the code
    in this class is optimized to be ~20X faster than the code in `Filtering`. This
    assumes we have the log model.

    It might seem unneccessary to have so many local attributes, but it speeds up
    the inference considerably if we index a value from an array once and store it
    as a float instead of repeatedly indexing the array - the difference in
    reality is super small but we do this so often that it adds up to ~40% speedup
    as supposed to not doing it - so we do this as often as possible - This speedup
    is even greater for indexing keys of dictionaries and getting parameters of objects.

    General efficiency speedups
    ---------------------------
    All of these are done relative to a non-optimized filtering implementation
        - Specialized sampling and logpdf functions. About 95% faster than
          scipy or numpy functions. All of these add up to a ~35% speed up
        - Explicit function definitions:
          instead of doing `self.prior.logpdf(...)`, we do
          `pl.random.normal.logpdf(...)`, about a ~10% overall speedup
        - Precomputation of values so that the least amout of computation
          is done on a data level - All of these add up to a ~25% speed up
    Benchmarked on a MacPro

    Non/trivial efficiency speedups w.r.t. non-multiprocessed filtering
    -------------------------------------------------------------------
    All of the speed ups are done relative to the current implementation
    in non-multiprocessed filtering.
        - Whenever possible we replace a 2D variable like `self.x` with a
          `curr_x`, which is 1D, because indexing a 1D array is 10-20% faster
          than a 2D array. All of these add up to a ~8% speed up
        - We only every compute the forward dynamics (not the reverse), because
          we can use the forward of the previous timepoint as the reverse for
          the next timepoint. This is about 45% faster and adds up to a ~40%
          speedup.
        - Whenever possible, we replace a dictionary like `read_depths`
          with a float because indexing a dict is 12-20% slower than a 1D
          array. All these add up to a ~7% speed up
        - We precompute AS MUCH AS POSSIBLE in `update` and in `initialize`,
          even simple this as `self.curr_tidx_minus_1`: all of these add up
          to about ~5% speedup
        - If an attribute of a class is being referenced more than once in
          a subroutine, we "get" it by making it a local variable. Example:
          `tidx = self.tidx`. This has about a 5% speed up PER ADDITIONAL
          CALL within the subroutine. All of these add up to ~2.5% speed up.
        - If an indexed value gets indexed more than once within a subroutine,
          we "get" the value by making it a local variable. All of these
          add up to ~4% speed up.
        - We "get" all of the means and stds of the qPCR data-structures so we
          do not reference an object. This is about 22% faster and adds up
          to a ~3% speed up.
    Benchmarked on a MacPro
    '''
    def __init__(self):
        '''Set all local variables to None
        '''
        return

    def initialize(self, times, qpcr_log_measurements, reads, there_are_intermediate_timepoints,
        there_are_perturbations, pv_global, x_prior_mean,
        x_prior_std, tune, delay, end_iter, proposal_init_scale, a0, a1, x, calculate_qpcr_loglik,
        pert_starts, pert_ends, ridx, subjname, h5py_xname, target_acceptance_rate,
        zero_inflation_transition_policy):
        '''Initialize the object at the beginning of the inference

        n_o = Number of ASVs
        n_gT = Number of given time points
        n_T = Total number of time points, including intermediate
        n_P = Number of Perturbations

        Parameters
        ----------
        times : np.array((n_T, ))
            Times for each of the time points
        qpcr_log_measurements : dict(t -> np.ndarray(float))
            These are the qPCR observations for every timepoint in log space.
        reads : dict (float -> np.ndarray((n_o, )))
            The counts for each of the given timepoints. Each value is an
            array for the counts for each of the ASVs
        there_are_intermediate_timepoints : bool
            If True, then there are intermediate timepoints, else there are only
            given timepoints
        there_are_perturbations : bool
            If True, that means there are perturbations, else there are no
            perturbations
        pv_global : bool
            If True, it means the process variance is is global for each ASV. If
            False it means that there is a separate `pv` for each ASV
        pv : float, np.ndarray
            This is the process variance value. This is a float if `pv_global` is True
            and it is an array if `pv_global` is False.
        x_prior_mean, x_prior_std : numeric
            This is the prior mean and std for `x` used when sampling the reverse for
            the first timepoint
        tune : int
            How often we should update the proposal for each ASV
        delay : int
            How many MCMC iterations we should delay the start of updating
        end_iter : int
            What iteration we should stop updating the proposal
        proposal_init_scale : float
            Scale to multiply the initial covariance of the poposal
        a0, a1 : floats
            These are the negative binomial dispersion parameters that specify how
            much noise there is in the counts
        x : np.ndarray((n_o, n_T))
            This is the x initialization
        pert_starts, pert_ends : np.ndarray((n_P, ))
            The starts and ends for each one of the perturbations
        ridx : int
            This is the replicate index that this object corresponds to
        subjname : str
            This is the name of the replicate
        h5py_xname : str
            This is the name for the x in the h5py object
        target_acceptance_rate : float
            This is the target acceptance rate for each point
        calculate_qpcr_loglik : bool
            If True, calculate the loglikelihood of the qPCR measurements during the proposal
        '''
        self.subjname = subjname
        self.h5py_xname = h5py_xname
        self.target_acceptance_rate = target_acceptance_rate
        self.zero_inflation_transition_policy = zero_inflation_transition_policy

        self.times = times
        self.qpcr_log_measurements = qpcr_log_measurements
        self.reads = reads
        self.there_are_intermediate_timepoints = there_are_intermediate_timepoints
        self.there_are_perturbations = there_are_perturbations
        self.pv_global = pv_global
        if not pv_global:
            raise TypeError('Filtering with MP not implemented for non global process variance')
        self.x_prior_mean = x_prior_mean
        self.x_prior_std = x_prior_std
        self.tune = tune
        self.delay = 0
        self.end_iter = end_iter
        self.proposal_init_scale = proposal_init_scale
        self.a0 = a0
        self.a1 = a1
        self.n_asvs = x.shape[0]
        self.n_timepoints = len(times)
        self.n_timepoints_minus_1 = len(times)-1
        self.logx = np.log(x)
        self.x = x

        # Get the perturbations for this subject
        self.pert_starts = []
        self.pert_ends = []
        for pidx in range(len(pert_starts)):
            self.pert_starts.append(pert_starts[pidx][subjname])
            self.pert_ends.append(pert_ends[pidx][subjname])

        self.total_n_points = self.x.shape[0] * self.x.shape[1]
        self.ridx = ridx
        self.calculate_qpcr_loglik = calculate_qpcr_loglik

        self.sample_iter = 0
        self.n_data_points = self.x.shape[0] * self.x.shape[1]

        # latent state
        self.sum_q = np.sum(self.x, axis=0)
        shape = (self.tune, ) + self.x.shape
        self.trace_iter = 0

        # proposal
        self.proposal_std = np.log(1.5) #np.log(3)
        self.acceptances = 0
        self.n_props_total = 0
        self.n_props_local = 0
        self.total_acceptances = 0
        self.add_trace = True

        # Intermediate timepoints
        if self.there_are_intermediate_timepoints:
            self.is_intermediate_timepoint = {}
            self.data_loglik = self.data_loglik_w_intermediates
            for t in self.times:
                self.is_intermediate_timepoint[t] = t not in self.reads
        else:
            self.data_loglik = self.data_loglik_wo_intermediates

        # Reads
        self.read_depths = {}
        for t in self.reads:
            self.read_depths[t] = float(np.sum(self.reads[t]))

        # t
        self.dts = np.zeros(self.n_timepoints_minus_1)
        self.sqrt_dts = np.zeros(self.n_timepoints_minus_1)
        for k in range(self.n_timepoints_minus_1):
            self.dts[k] = self.times[k+1] - self.times[k]
            self.sqrt_dts[k] = np.sqrt(self.dts[k])
        self.t2tidx = {}
        for tidx, t in enumerate(self.times):
            self.t2tidx[t] = tidx

        self.cnt_accepted_times = np.zeros(len(self.times))

        # Perturbations
        # -------------
        # in_pert_transition : np.ndarray(dtype=bool)
        #   This is a bool array where if it is a true it means that the
        #   forward and reverse growth rates are different
        # fully_in_pert : np.ndarray(dtype=int)
        #   This is an int-array where it tells you which perturbation you are fully in
        #   (the forward and reverse growth rates are the same but not the default).
        #   If there is no perturbation then the value is -1. If it is not -1, then the
        #   number corresponds to what perturbation index you are in.
        #
        # Edge cases
        # ----------
        #   * missing data for start
        #       There could be a situation where there was no sample collection
        #       on the day that they started a perturbation. In this case we
        #       assume that the next time point is the `start` of the perturbation.
        #       i.e. the next time point is the perturbation transition.
        #   * missing data for end
        #       There could be a situation where no sample was collected when the
        #       perturbation ended. In this case we assume that the pervious time
        #       point was the end of the perturbation.
        if self.there_are_perturbations:
            self.in_pert_transition = np.zeros(self.n_timepoints, dtype=bool)
            self.fully_in_pert = np.ones(self.n_timepoints, dtype=int) * -1
            for pidx, t in enumerate(self.pert_starts):
                if t == self.times[-1] or t == self.times[0]:
                    raise ValueError('The code right now does not support either a perturbation that ' \
                        'started on the first day or ended on the last day. The code where this is ' \
                        'incompatible is when we checking if we are in a perturbation transition')
                if t > np.max(self.times):
                    continue
                if t not in self.t2tidx:
                    # Use the next time point
                    tidx = np.searchsorted(self.times, t)
                else:
                    tidx = self.t2tidx[t]
                self.in_pert_transition[tidx] = True
            for pidx, t in enumerate(self.pert_ends):
                if t == self.times[-1] or t == self.times[0]:
                    raise ValueError('The code right now does not support either a perturbation that ' \
                        'started on the first day or ended on the last day. The code where this is ' \
                        'incompatible is when we checking if we are in a perturbation transition')
                if t < np.min(self.times) or t > np.max(self.times):
                    continue
                if t not in self.t2tidx:
                    # Use the previous time point
                    tidx = np.searchsorted(self.times, t) - 1
                else:
                    tidx = self.t2tidx[t]
                self.in_pert_transition[tidx] = True

            # check if anything is weird
            if np.sum(self.in_pert_transition) % 2 != 0:
                raise ValueError('The number of in_pert_transition periods must be even ({})' \
                    '. There is either something wrong with the data (start and end day are ' \
                    'the same) or with the algorithm ({})'.format(
                        np.sum(self.in_pert_transition),
                        self.in_pert_transition))

            # Make the fully in perturbation times
            for pidx in range(len(self.pert_ends)):
                try:
                    start_tidx = self.t2tidx[self.pert_starts[pidx]] + 1
                    end_tidx = self.t2tidx[self.pert_ends[pidx]]
                except:
                    # This means there is a missing datapoint at either the
                    # start or end of the perturbation
                    start_t = self.pert_starts[pidx]
                    end_t = self.pert_ends[pidx]
                    start_tidx = np.searchsorted(self.times, start_t)
                    end_tidx = np.searchsorted(self.times, end_t) - 1

                self.fully_in_pert[start_tidx:end_tidx] = pidx
    
    # @profile
    def persistent_run(self, growth, self_interactions, pv, interactions,
        perturbations, qpcr_variances, zero_inflation_data):
        '''Run an update of the values for a single gibbs step for all of the data points
        in this replicate

        Parameters
        ----------
        growth : np.ndarray((n_asvs, ))
            Growth rates for each ASV
        self_interactions : np.ndarray((n_asvs, ))
            Self-interactions for each ASV
        pv : numeric, np.ndarray
            This is the process variance
        interactions : np.ndarray((n_asvs, n_asvs))
            These are the ASV-ASV interactions
        perturbations : np.ndarray((n_perturbations, n_asvs))
            Perturbation values in the right perturbation order, per ASV
        zero_inflation : np.ndarray
            These are the points that are delibertly pushed down to zero
        qpcr_variances : np.ndarray
            These are the sampled qPCR variances as an array - they are in
            time order

        Returns
        -------
        (int, np.ndarray, float)
            1 This is the replicate index
            2 This is the updated latent state for logx
            3 This is the acceptance rate for this past update.
        '''
        self.master_growth_rate = growth

        if self.sample_iter < self.delay:
            self.sample_iter += 1
            return self.ridx, self.x, np.nan

        self.update_proposals()
        self.n_accepted_iter = 0
        self.pv = pv
        self.pv_std = np.sqrt(pv)
        self.qpcr_stds = np.sqrt(qpcr_variances[self.ridx])
        self.qpcr_stds_d = {}
        # self.zero_inflation_data = zero_inflation_data[self.ridx]
        self.zero_inflation_data = None

        for tidx,t in enumerate(self.qpcr_log_measurements):
            self.qpcr_stds_d[t] = self.qpcr_stds[tidx]

        if self.there_are_perturbations:
            self.growth_rate_non_pert = growth.ravel()
            self.growth_rate_on_pert = growth.reshape(-1,1) * (1 + perturbations)
                
        # Go through each randomly ASV and go in time order
        oidxs = npr.permutation(self.n_asvs)
        # print('===============================')
        # print('===============================')
        # print('ridx', self.ridx)
        for oidx in oidxs:

            # Set the necessary global parameters
            self.oidx = oidx
            self.curr_x = self.x[oidx, :]
            self.curr_logx = self.logx[oidx, :]
            self.curr_interactions = interactions[oidx, :]
            self.curr_self_interaction = self_interactions[oidx]
            # self.curr_zero_inflation = self.zero_inflation[oidx, :]

            if self.pv_global:
                self.curr_pv_std = self.pv_std
            else:
                self.curr_pv_std = self.pv_std[oidx]

            # Set for first time point
            self.tidx = 0
            self.set_attrs_for_timepoint()
            self.forward_loglik = self.default_forward_loglik
            self.reverse_loglik = self.first_timepoint_reverse
            # Calculate A matrix for forward
            self.forward_interaction_vals = np.nansum(self.x[:, self.tidx] * self.curr_interactions)
            self.update_single()
            self.reverse_loglik = self.default_reverse_loglik
            # Set for middle timepoints
            for tidx in range(1, self.n_timepoints-1):
                # Check if it needs to be zero inflated
                # if not self.curr_zero_inflation[tidx]:
                #     raise NotImplementedError('Zero inflation not implemented for logmodel')

                self.tidx = tidx
                self.set_attrs_for_timepoint()

                # Calculate A matrix for forward and reverse
                # Set the reverse of the current time step to the forward of the previous
                self.reverse_interaction_vals = self.forward_interaction_vals #np.sum(self.x[:, self.prev_tidx] * self.curr_interactions)
                self.forward_interaction_vals = np.nansum(self.x[:, self.tidx] * self.curr_interactions)

                # Run single update
                self.update_single()

            # Set for last timepoint
            self.tidx = self.n_timepoints_minus_1
            self.set_attrs_for_timepoint()
            self.forward_loglik = self.last_timepoint_forward
            # Calculate A matrix for reverse
            # Set the reverse of the current time step to the forward of the previous
            self.reverse_interaction_vals = self.forward_interaction_vals # np.sum(self.x[:, self.prev_tidx] * self.curr_interactions)
            self.update_single()

            # if self.sample_iter == 4:
            # sys.exit()

        self.sample_iter += 1
        if self.add_trace:
            self.trace_iter += 1

        # print(self.cnt_accepted_times/self.sample_iter)

        return self.ridx, self.x, self.n_accepted_iter/self.n_data_points

    def set_attrs_for_timepoint(self):
        self.prev_tidx = self.tidx-1
        self.next_tidx = self.tidx+1
        self.forward_growth_rate = self.master_growth_rate[self.oidx]
        self.reverse_growth_rate = self.master_growth_rate[self.oidx]

        if self.there_are_intermediate_timepoints:
            if not self.is_intermediate_timepoint[self.times[self.tidx]]:
                # It is not intermediate timepoints - we need to get the data
                t = self.times[self.tidx]
                self.curr_reads = self.reads[t][self.oidx]
                self.curr_read_depth = self.read_depths[t]
                self.curr_qpcr_log_measurements = self.qpcr_log_measurements[t]
                self.curr_qpcr_std = self.qpcr_stds_d[t]
        else:
            t = self.times[self.tidx]
            self.curr_reads = self.reads[t][self.oidx]
            self.curr_read_depth = self.read_depths[t]
            self.curr_qpcr_log_measurements = self.qpcr_log_measurements[t]
            self.curr_qpcr_std = self.qpcr_stds_d[t]

        # Set perturbation growth rates
        if self.there_are_perturbations:
            if self.in_pert_transition[self.tidx]:
                if self.fully_in_pert[self.tidx-1] != -1:
                    # If the previous time point is in the perturbation, that means
                    # we are going out of the perturbation
                    # self.forward_growth_rate = self.master_growth_rate[self.oidx]
                    pidx = self.fully_in_pert[self.tidx-1]
                    self.reverse_growth_rate = self.growth_rate_on_pert[self.oidx,pidx]
                else:
                    # Else we are going into a perturbation
                    # self.reverse_growth_rate = self.master_growth_rate[self.oidx]
                    pidx = self.fully_in_pert[self.tidx+1]
                    self.forward_growth_rate = self.growth_rate_on_pert[self.oidx,pidx]
            elif self.fully_in_pert[self.tidx] != -1:
                pidx = self.fully_in_pert[self.tidx]
                self.forward_growth_rate = self.growth_rate_on_pert[self.oidx,pidx]
                self.reverse_growth_rate = self.forward_growth_rate

    # @profile
    def update_single(self):
        '''Update a single oidx, tidx
        '''
        tidx = self.tidx
        oidx = self.oidx

        # Check if we should update the zero inflation policy
        if self.zero_inflation_transition_policy is not None:
            if self.zero_inflation_transition_policy == 'ignore':
                if not self.zero_inflation_data[oidx,tidx]:
                    self.x[oidx, tidx] = np.nan
                    self.logx[oidx, tidx] = np.nan
                    return
                else:
                    if tidx < self.zero_inflation_data.shape[1]-1:
                        do_forward = self.zero_inflation_data[oidx, tidx+1]
                    else:
                        do_forward = True
                    if tidx > 0:
                        do_reverse = self.zero_inflation_data[oidx, tidx-1]
                    else:
                        do_reverse = True
            else:
                raise NotImplementedError('Not Implemented')
        else:
            do_forward = True
            do_reverse = True

        # t = self.times[self.tidx]
        # # proposal
        # mu1 = self.curr_logx[tidx]
        # rel = self.reads[t][oidx]/self.read_depths[t]
        # if rel == 0:
        #     rel = 1e-5
        # mu2 = np.log(rel*np.exp(self.curr_qpcr_loc + (self.curr_qpcr_scale/2)))

        # var1 = self.proposal_std[(tidx, oidx)]**2
        # var2 = (self.curr_qpcr_scale)**2
        # mu,var = prod_gaussians(means=[mu1,mu2], variances=[var1,var2])

        try:
            logx_new = pl.random.misc.fast_sample_normal(
                self.curr_logx[tidx],
                self.proposal_std)
        except:
            print('mu', self.curr_logx[tidx])
            print('std', self.proposal_std)
            raise
        # try:
        #     logx_new = pl.random.misc.fast_sample_normal(
        #         mu, np.sqrt(var))
        # except:
        #     print('mu', mu)
        #     print('std', np.sqrt(var))
        #     raise

        x_new = np.exp(logx_new)
        prev_logx_value = self.curr_logx[tidx]
        prev_x_value = self.curr_x[tidx]

        # print('prex_x', prev_x_value, np.exp(prev_logx_value))
        # print('prev_logx', prev_logx_value)

        # if tidx == 5:
        #     print('\ntidx', tidx)
        #     print('oidx', oidx)
        #     print('t', self.times[self.tidx])
        #     print('curr_logx', self.curr_logx[tidx])
        #     # print('curr_logx', self.curr_logx[tidx])
        #     # print('mu1, mu2', mu1, mu2)
        #     # print('mu', mu)
        #     # print('var1, var2', var1,var2)
        #     # print('var',var)
        #     print('prop_logx', logx_new)
        #     # print('start perts', self.pert_starts)
        #     # print('end perts', self.pert_ends)
        #     # print('in perturbation transition?', self.in_pert_transition[tidx])
        #     # print('fully in pert?', self.fully_in_pert[self.tidx])
        #     print('forward growth', self.forward_growth_rate)

        if do_forward:
            prev_aaa = self.forward_loglik()
        else:
            prev_aaa = 0
        if do_reverse:
            prev_bbb = self.reverse_loglik()
        else:
            prev_bbb = 0
        prev_ddd = self.data_loglik()

        # if tidx == 5:
        #     print('\nold')
        #     print('forward ll', aaa)
        #     print('reverse ll', bbb)
        #     print('data ll', ddd)

        l_old = prev_aaa + prev_bbb + prev_ddd

        self.curr_x[tidx] = x_new
        self.curr_logx[tidx] = logx_new
        self.sum_q[tidx] = self.sum_q[tidx] - prev_x_value + x_new

        if do_forward:
            new_aaa = self.forward_loglik()
        else:
            new_aaa = 0
        if do_reverse:
            new_bbb = self.reverse_loglik()
        else:
            new_bbb = 0
        new_ddd = self.data_loglik()

        # if tidx == 5:
        #     print('\nnew')
        #     print('forward ll', aaa)
        #     print('reverse ll', bbb)
        #     print('data ll', ddd)
        #     print('\nold x value', prev_x_value)
        #     print('old logx value', prev_logx_value)
        #     print('proposal std', self.proposal_std[(tidx, oidx)])
        #     print('new x value', x_new)
        #     print('new logx value', logx_new)

        l_new = new_aaa + new_bbb + new_ddd
        r_accept = l_new - l_old

        # if tidx == 0:
        #     print('\n\noidx {} diff lls:'.format(oidx), r_accept)
        #     print('\tforward', new_aaa - prev_aaa)
        #     print('\treverse', new_bbb - prev_bbb)
        #     print('\tdata', new_ddd - prev_ddd)

        # if tidx == 5:
        #     print('r_accept', r_accept)
        r = pl.random.misc.fast_sample_standard_uniform()
        if math.log(r) > r_accept:
            # print('reject')
            # reject
            self.sum_q[tidx] = self.sum_q[tidx] + prev_x_value - x_new
            self.curr_x[tidx] = prev_x_value
            self.curr_logx[tidx] = prev_logx_value
        else:
            # print('accept')
            self.x[oidx, tidx] = x_new
            self.logx[oidx, tidx] = logx_new
            self.acceptances += 1
            self.total_acceptances += 1
            self.n_accepted_iter += 1

        self.n_props_local += 1
        self.n_props_total += 1

    def update_proposals(self):
        '''Update the proposal if necessary
        '''
        if self.sample_iter > self.end_iter:
            self.add_trace = False
            return
        if self.sample_iter == 0:
            return
        if self.trace_iter - self.delay == self.tune  and self.sample_iter - self.delay > 0:

            # Adjust
            acc_rate = self.acceptances/self.n_props_total
            if acc_rate < 0.1:
                logging.debug('Very low acceptance rate, scaling down past covariance')
                self.proposal_std *= 0.01
            elif acc_rate < self.target_acceptance_rate:
                self.proposal_std /= np.sqrt(1.5)
            else:
                self.proposal_std *= np.sqrt(1.5)
            
            self.acceptances = 0
            self.n_props_local = 0

    def last_timepoint_forward(self):
        return 0

    # @profile
    def default_forward_loglik(self):
        '''From the current timepoint (tidx) to the next timepoint (tidx+1)
        '''
        logmu = self.compute_dynamics(
            tidx=self.tidx,
            Axj=self.forward_interaction_vals,
            a1=self.forward_growth_rate)

        try:
            return pl.random.normal.logpdf(
                value=self.curr_logx[self.next_tidx], 
                mean=logmu, std=self.curr_pv_std*self.sqrt_dts[self.tidx])
        except:
            return _normal_logpdf(
                value=self.curr_logx[self.next_tidx], 
                mean=logmu, std=self.curr_pv_std*self.sqrt_dts[self.tidx])

    def first_timepoint_reverse(self):
        # sample from the prior
        try:
            return pl.random.normal.logpdf(value=self.curr_logx[self.tidx],
                mean=self.x_prior_mean, std=self.x_prior_std)
        except:
            return _normal_logpdf(value=self.curr_logx[self.tidx],
                mean=self.x_prior_mean, std=self.x_prior_std)

    # @profile
    def default_reverse_loglik(self):
        '''From the previous timepoint (tidx-1) to the current time point (tidx)
        '''
        logmu = self.compute_dynamics(
            tidx=self.prev_tidx,
            Axj=self.reverse_interaction_vals,
            a1=self.reverse_growth_rate)

        try:
            return pl.random.normal.logpdf(value=self.curr_logx[self.tidx], 
                mean=logmu, std=self.curr_pv_std*self.sqrt_dts[self.prev_tidx])
        except:
            return _normal_logpdf(value=self.curr_logx[self.tidx], 
                mean=logmu, std=self.curr_pv_std*self.sqrt_dts[self.prev_tidx])

    def data_loglik_w_intermediates(self):
        '''data loglikelihood w/ intermediate timepoints
        '''
        if self.is_intermediate_timepoint[self.times[self.tidx]]:
            return 0
        else:
            return self.data_loglik_wo_intermediates()

    # @profile
    def data_loglik_wo_intermediates(self):
        '''data loglikelihood with intermediate timepoints
        '''
        sum_q = self.sum_q[self.tidx]
        log_sum_q = math.log(sum_q)
        rel = self.curr_x[self.tidx] / sum_q

        try:
            negbin = negbin_loglikelihood_MH_condensed(
                k=self.curr_reads,
                m=self.curr_read_depth * rel,
                dispersion=self.a0/rel + self.a1)
        except:
            negbin = negbin_loglikelihood_MH_condensed_not_fast(
                k=self.curr_reads,
                m=self.curr_read_depth * rel,
                dispersion=self.a0/rel + self.a1)

        qpcr = 0
        if self.calculate_qpcr_loglik:
            for qpcr_val in self.curr_qpcr_log_measurements:
                a = pl.random.normal.logpdf(value=qpcr_val, mean=log_sum_q, std=self.curr_qpcr_std)
                qpcr += a

        # tidx = self.tidx
        # if True: #tidx in [3,6,7]:
        #     print('\n\nData, tidx', tidx)
        #     print('sum_q:', sum_q)
        #     print('rel:', rel)
        #     print('qpcr: {}\n\tvalue: {}\n\tmean: {}\n\tstd: {}'.format(
        #         qpcr, self.curr_qpcr_loc,
        #         sum_q,
        #         self.curr_qpcr_scale))
        #     print('neg_bin: {}\n\tk: {}\n\tm: {}\n\tdispersion: {}'.format( 
        #         negbin, self.curr_reads, self.curr_read_depth * rel,
        #         self.a0/rel + self.a1))
            
        #     print('data\n\tcurr_x: {}, {}\n\tcurr_logx: {}'.format(
        #         self.curr_x[self.tidx], 
        #         np.exp(self.curr_logx[self.tidx]),
        #         self.curr_logx[self.tidx]))
        return negbin + qpcr

    def compute_dynamics(self, tidx, Axj, a1):
        '''Compute dynamics going into tidx+1

        a1 : growth rates and perturbations (if necessary)
        Axj : cluster interactions with the other abundances already multiplied
        tidx : time index

        Zero-inflation
        --------------
        When we get here, the current asv at `tidx` is not a structural zero, but 
        there might be other bugs in the system that do have a structural zero there.
        Thus we do nan adds
        '''
        logxi = self.curr_logx[tidx]
        xi = self.curr_x[tidx]

        # print('dynamics')
        # print('xi*a1', xi*a1* self.dts[tidx])
        # print('xi*xi*self.curr_self_interaction', xi*xi*self.curr_self_interaction* self.dts[tidx])
        # print('xi*Axj', xi*Axj* self.dts[tidx])

        # compute dynamics
        return logxi + (a1 - xi*self.curr_self_interaction + Axj) * self.dts[tidx]


class ZeroInflation(pl.graph.Node):
    '''This is the posterior distribution for the zero inflation model. These are used
    to learn when the model should use use the data and when it should not. We do not need
    to trace this object because we set the structural zeros to nans in the trace for 
    filtering.

    TODO: Parallel version of the class
    '''

    def __init__(self, **kwargs):
        '''
        Parameters
        ----------
        mp : str
            This is the type of parallelization to use. This is not implemented yet.
        '''
        kwargs['name'] = STRNAMES.ZERO_INFLATION
        pl.graph.Node.__init__(self, **kwargs)
        self.value = []
        self._strr = 'NA'

        for ridx in range(self.G.data.n_replicates):
            n_timepoints = self.G.data.n_timepoints_for_replicate[ridx]
            self.value.append(np.ones(shape=(len(self.G.data.asvs), n_timepoints), dtype=bool))

    def reset_value_size(self):
        '''Change the size of the trajectory when we set the intermediate timepoints
        '''
        n_asvs = self.G.data.n_asvs
        for ridx in range(len(self.value)):
            n_timepoints = self.G.data.n_timepoints_for_replicate[ridx]
            self.value[ridx] = np.ones((n_asvs, n_timepoints), dtype=bool)

    def __str__(self):
        return self._strr

    def initialize(self, value_option, delay=0):
        '''Initialize the values. Right now this is static and we are not learning this so
        do not do anything fancy

        Parameters
        ----------

        delay : None, int
            How much to delay starting the sampling
        '''
        if delay is None:
            delay = 0
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        self.delay = delay

        if value_option in [None, 'auto']:
            # Set everything to on
            self.value = []
            for ridx in range(self.G.data.n_replicates):
                n_timepoints = self.G.data.n_timepoints_for_replicate[ridx]
                self.value.append(np.ones(
                    shape=(len(self.G.data.asvs), n_timepoints), dtype=bool))
            turn_on = None
            turn_off = None

        elif value_option == 'mdsine-cdiff':
            # Set everything to on except for cdiff before day 28 for every subject
            self.value = []
            for ridx in range(self.G.data.n_replicates):
                n_timepoints = self.G.data.n_timepoints_for_replicate[ridx]
                self.value.append(np.ones(
                    shape=(len(self.G.data.asvs), n_timepoints), dtype=bool))

            # Get cdiff
            cdiff_idx = self.G.data.asvs['Clostridium-difficile'].idx
            turn_off = []
            turn_on = []
            for ridx in range(self.G.data.n_replicates):
                for tidx, t in enumerate(self.G.data.times[ridx]):
                    for oidx in range(len(self.G.data.asvs)):
                        if t < 28 and oidx == cdiff_idx:
                            self.value[ridx][cdiff_idx, tidx] = False
                            turn_off.append((ridx, tidx, cdiff_idx))
                        else:
                            turn_on.append((ridx, tidx, oidx))

        else:
            raise ValueError('`value_option` ({}) not recognized'.format(value_option))

        self.G.data.set_zero_inflation(turn_on=turn_on, turn_off=turn_off)
                

# Interactions
# ------------
class PriorVarInteractions(pl.variables.SICS):
    '''This is the posterior of the prior variance of regression coefficients
    for the interaction (off diagonal) variables
    '''
    def __init__(self, prior, value=None, **kwargs):

        kwargs['name'] = STRNAMES.PRIOR_VAR_INTERACTIONS
        pl.variables.SICS.__init__(self, value=value,
            dtype=float, **kwargs)
        self.add_prior(prior)

    def initialize(self, value_option, dof_option, scale_option, value=None,
        mean_scaling_factor=None, dof=None, scale=None, delay=0):
        '''Initialize the hyperparameters of the self interaction variance based on the
        passed in option

        Parameters
        ----------
        value_option : str
            - Initialize the value based on the specified option
            - Options
                - 'manual'
                    - Set the value manually, `value` must also be specified
                - 'auto', 'prior-mean'
                    - Set the value to the mean of the prior
        scale_option : str
            - Initialize the scale of the prior
            - Options
                - 'manual'
                    - Set the value manually, `scale` must also be specified
                - 'auto', 'same-as-aii'
                    - Set the mean the same as the self-interactions
        dof_option : str
            Initialize the dof of the parameter
            Options:
                'manual': Set the value with the parameter `dof`
                'diffuse': Set the value to 2.01
                'strong': Set the valuye to the expected number of interactions
                'auto': Set to diffuse
        dof, scale : int, float
            - User specified values
            - Only necessary if `hyperparam_option` == 'manual'
        '''
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        self.delay = delay

        self.interactions = self.G[REPRNAMES.CLUSTER_INTERACTION_VALUE]

        if not pl.isstr(dof_option):
            raise TypeError('`dof_option` ({}) must be a str'.format(type(dof_option)))
        if dof_option == 'manual':
            if not pl.isnumeric(dof):
                raise TypeError('`dof` ({}) must be a numeric'.format(type(dof)))
            if dof < 0:
                raise ValueError('`dof` ({}) must be > 0 for it to be a valid prior'.format(dof))
        elif dof_option in ['diffuse', 'auto']:
            dof = 2.01
        elif dof_option == 'strong':
            N = expected_n_clusters(G=self.G)
            dof = N * (N - 1)
        else:
            raise ValueError('`dof_option` ({}) not recognized'.format(dof_option))
        self.prior.dof.override_value(dof)

        if not pl.isstr(scale_option):
            raise TypeError('`scale_option` ({}) must be a str'.format(type(scale_option)))
        if scale_option == 'manual':
            if not pl.isnumeric(scale):
                raise TypeError('`scale` ({}) must be a numeric'.format(type(scale)))
            if scale < 0:
                raise ValueError('`scale` ({}) must be > 0 for it to be a valid prior'.format(scale))
        elif scale_option in ['auto', 'same-as-aii']:
            mean = self.G[REPRNAMES.PRIOR_VAR_SELF_INTERACTIONS].prior.mean()
            scale = mean * (self.prior.dof.value - 2) /(self.prior.dof.value)
        else:
            raise ValueError('`scale_option` ({}) not recognized'.format(scale_option))
        self.prior.scale.override_value(scale)

        if not pl.isstr(value_option):
            raise TypeError('`value_option` ({}) must be a str'.format(type(value_option)))
        if value_option == 'manual':
            if not pl.isnumeric(value):
                raise ValueError('`value` ({}) must be numeric (float,int)'.format(value.__class__))
            self.value = value
        elif value_option in ['auto', 'prior-mean']:
            if not pl.isnumeric(mean_scaling_factor):
                raise ValueError('`mean_scaling_factor` ({}) must be a numeric type ' \
                    '(float,int)'.format(mean_scaling_factor.__class__))
            self.value = self.prior.mean()
        else:
            raise ValueError('`value_option` ({}) not recognized'.format(value_option))

        logging.info('Prior Variance Interactions initialization results:\n' \
            '\tprior dof: {}\n' \
            '\tprior scale: {}\n' \
            '\tvalue: {}'.format(
                self.prior.dof.value, self.prior.scale.value, self.value))

    # @profile
    def update(self):
        '''Calculate the posterior of the prior variance
        '''
        if self.sample_iter < self.delay:
            return

        x = self.interactions.obj.get_values(use_indicators=True)
        mu = self.G[REPRNAMES.PRIOR_MEAN_INTERACTIONS].value

        se = np.sum(np.square(x - mu))
        n = len(x)

        self.dof.value = self.prior.dof.value + n
        self.scale.value = ((self.prior.scale.value * self.prior.dof.value) + \
           se)/self.dof.value
        self.sample()

    def visualize(self, path, f, section='posterior'):
        return _scalar_visualize(self, path=path, f=f, section=section)


class PriorMeanInteractions(pl.variables.Normal):
    '''This is the posterior mean for the interactions
    '''

    def __init__(self, prior, **kwargs):
        kwargs['name'] = STRNAMES.PRIOR_MEAN_INTERACTIONS
        pl.variables.Normal.__init__(self, mean=None, var=None, dtype=float, **kwargs)
        self.add_prior(prior)

    def __str__(self):
        # If this fails, it is because we are dividing by 0 sampler_iter
        # If which case we just return the value 
        try:
            s = 'Value: {}, Acceptance rate: {}'.format(
                self.value, np.mean(self.acceptances[
                    np.max([self.sample_iter-50, 0]):self.sample_iter]))
        except:
            s = str(self.value)
        return s

    def initialize(self, value_option, mean_option, var_option, value=None, 
        mean=None, var=None, delay=0):
        '''Initialize the hyperparameters

        Parameters
        ----------
        value_option : str
            How to set the value. Options:
                'zero'
                    Set to zero
                'prior-mean', 'auto'
                    Set to the mean of the prior
                'manual'
                    Specify with the `value` parameter
        mean_option : str
            How to set the mean of the prior
                'zero', 'auto'
                    Set to zero
                'manual'
                    Set with the `mean` parameter
        var_option : str
            'same-as-aii', 'auto'
                Set as the same variance as the self-interactions
            'manual'
                Set with the `var` parameter
        value, mean, var : float
            These are only necessary if we specify manual for any of the other 
            options
        delay : int
            How much to delay the start of the update during inference
        '''
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        self.delay = delay

        # Set the mean
        if not pl.isstr(mean_option):
            raise TypeError('`mean_option` ({}) must be a str'.format(type(mean_option)))
        if mean_option == 'manual':
            if not pl.isnumeric(mean):
                raise TypeError('`mean` ({}) must be a numeric'.format(type(mean)))
        elif mean_option in ['zero', 'auto']:
            mean = 0
        else:
            raise ValueError('`mean_option` ({}) not recognized'.format(mean_option))
        self.prior.mean.override_value(mean)

        # Set the variance
        if not pl.isstr(var_option):
            raise TypeError('`var_option` ({}) must be a str'.format(type(var_option)))
        if var_option == 'manual':
            if not pl.isnumeric(var):
                raise TypeError('`var` ({}) must be a numeric'.format(type(var)))
            if var <= 0:
                raise ValueError('`var` ({}) must be positive'.format(var))
        elif var_option in ['same-as-aii', 'auto']:
            var = self.G[STRNAMES.PRIOR_VAR_SELF_INTERACTIONS].value
        else:
            raise ValueError('`var_option` ({}) not recognized'.format(var_option))
        self.prior.var.override_value(var)

        # Set the value
        if not pl.isstr(value_option):
            raise TypeError('`value_option` ({}) must be a str'.format(type(value_option)))
        if value_option == 'manual':
            if not pl.isnumeric(value):
                raise TypeError('`value` ({}) must be a numeric'.format(type(value)))
        elif value_option in ['prior-mean', 'auto']:
            value = self.prior.mean.value
        elif value_option == 'zero':
            value = 0
        else:
            raise ValueError('`value_option` ({}) not recognized'.format(value_option))
        self.value = value

    def update(self):
        '''Update using Gibbs sampling
        '''
        if self.sample_iter < self.delay:
            return

        if self.G[REPRNAMES.CLUSTER_INTERACTION_INDICATOR].num_pos_indicators == 0:
            # sample from the prior
            self.value = self.prior.sample()
            return

        x = self.G[REPRNAMES.CLUSTER_INTERACTION_VALUE].value
        prec = 1/self.G[REPRNAMES.PRIOR_VAR_INTERACTIONS].value

        prior_prec = 1/self.prior.var.value
        prior_mean = self.prior.mean.value

        self.var.value = 1/(prior_prec + (len(x)*prec))
        self.mean.value = self.var.value * ((prior_mean * prior_prec) + (np.sum(x)*prec))
        self.sample()

    def visualize(self, path, f, section='posterior'):
        return _scalar_visualize(self, path=path, f=f, section=section)


class ClusterInteractionValue(pl.variables.MVN):
    '''Interactions of Lotka-Voltera

    Since we initialize the interactions object in the `initialize` function,
    make sure that you have initialized the prior of the values of the interactions
    and of the indicators of the interactions before you call the initialization of
    this class
    '''
    def __init__(self, prior, clustering, **kwargs):
        kwargs['name'] = STRNAMES.CLUSTER_INTERACTION_VALUE
        pl.variables.MVN.__init__(self, dtype=float, **kwargs)
        self.set_value_shape(shape=(len(self.G.data.asvs),len(self.G.data.asvs)))
        self.add_prior(prior)
        self.clustering = clustering
        self.obj = pl.contrib.Interactions(
            clustering=self.clustering,
            use_indicators=True,
            name=STRNAMES.INTERACTIONS_OBJ, G=self.G,
            signal_when_clusters_change=False)
        self._strr = 'None'

    def __str__(self):
        return self._strr

    def __len__(self):
        # Return the number of on interactions
        return self.obj.num_pos_indicators()

    def set_values(self, *args, **kwargs):
        '''Set the values from an array
        '''
        self.obj.set_values(*args, **kwargs)

    def initialize(self, value_option, hyperparam_option=None, value=None,
        indicators=None, delay=0):
        '''Initialize the interactions object.

        Parameters
        ----------
        value_option : str
            This is how to initialize the values
            Options:
                'manual'
                    Set the values of the interactions manually. `value` and `indicators`
                    must also be specified. We assume the values are only set for when
                    `indicators` is True, and that the order of the `indicators` and `values`
                    correspond to how we iterate over the interactions
                    Example
                        3 Clusters
                        indicators = [True, False, False, True, False, True]
                        value = [0.2, 0.8, -0.35]
                'all-off', 'auto'
                    Set all of the interactions and the indicators to 0
                'all-on'
                    Set all of the indicators to on and all the values to 0
        delay : int
            How many MCMC iterations to delay starting to update
        See also
        --------
        `pylab.cluster.Interactions.set_values`
        '''
        self.obj.set_signal_when_clusters_change(True)
        self.G[REPRNAMES.INTERACTIONS_OBJ].value_initializer = self.prior.sample
        self.G[REPRNAMES.INTERACTIONS_OBJ].indicator_initializer = self.G[REPRNAMES.INDICATOR_PROB].prior.sample

        self._there_are_perturbations = self.G.perturbations is not None

        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        self.delay = delay

        if not pl.isstr(value_option):
            raise TypeError('`value_option` ({}) must be a str'.format(type(value_option)))
        if value_option in ['auto', 'all-off']:
            for interaction in self.obj:
                interaction.value = 0
                interaction.indicator = False
        elif value_option == 'all-on':
            for interaction in self.obj:
                interaction.value = 0
                interaction.indicator = True
        elif value_option == 'manual':
            if not np.all(pl.itercheck([value, indicators], pl.isarray)):
                raise TypeError('`value` ({}) and `indicators` ({}) must be arrays'.format(
                    type(value), type(indicators)))
            if len(value) != np.sum(indicators):
                raise ValueError('Length of `value` ({}) must equal the number of positive ' \
                    'values in `indicators` ({})'.format(len(value), np.sum(indicators)))
            if len(indicators) != self.obj.size:
                raise ValueError('The length of `indicators` ({}) must be the same as the ' \
                    'number of possible interactions ({})'.format(len(indicators), self.obj.size))
            ii = 0
            for i,interaction in enumerate(self.obj):
                interaction.indicator = indicators[i]
                if interaction.indicator:
                    interaction.value = value[ii]
                    ii += 1
        else:
            raise ValueError('`value_option` ({}) not recognized'.format(value_option))

        self._strr = str(self.obj.get_values(use_indicators=True))
        self.value = self.obj.get_values(use_indicators=True)

    def update(self):
        '''Update the values (where the indicators are positive) using a multivariate normal
        distribution - call this from regress coeff if you want to update the interactions
        conditional on all the other parameters.
        '''
        if self.obj.sample_iter < self.delay:
            return
        if self.obj.num_pos_indicators() == 0:
            # logging.info('No positive indicators, skipping')
            self._strr = '[]'
            return

        rhs = [
            REPRNAMES.CLUSTER_INTERACTION_VALUE]
        lhs = [
            REPRNAMES.GROWTH_VALUE,
            REPRNAMES.SELF_INTERACTION_VALUE]
        X = self.G.data.construct_rhs(keys=rhs)
        y = self.G.data.construct_lhs(keys=lhs,
            kwargs_dict={REPRNAMES.GROWTH_VALUE:{
                'with_perturbations':self._there_are_perturbations}})
        process_prec = self.G[REPRNAMES.PROCESSVAR].build_matrix(
            cov=False, sparse=True)
        prior_prec = build_prior_covariance(G=self.G, cov=False,
            order=rhs, sparse=True)

        pm = prior_prec @ (self.prior.mean.value * np.ones(prior_prec.shape[0]).reshape(-1,1))

        prec = X.T @ process_prec @ X + prior_prec
        cov = pinv(prec, self)
        mean = (cov @ (X.T @ process_prec.dot(y) + pm)).ravel()

        # print(np.hstack((y, self.G.data.lhs.vector.reshape(-1,1))))
        # for perturbation in self.G.perturbations:
        #     print()
        #     print(perturbation.magnitude.cluster_array())
        #     print(perturbation.indicator.cluster_array())

        self.mean.value = mean
        self.cov.value = cov
        value = self.sample()
        self.obj.set_values(arr=value, use_indicators=True)
        self.update_str()

        if np.any(np.isnan(self.value)):
            logging.critical('mean: {}'.format(self.mean.value))
            logging.critical('nan in cov: {}'.format(np.any(np.isnan(self.cov.value))))
            logging.critical('value: {}'.format(self.value))
            raise ValueError('`Values in {} are nan: {}'.format(self.name, self.value))

    def update_str(self):
        self._strr = str(self.obj.get_values(use_indicators=True))

    def set_trace(self):
        self.obj.set_trace()

    def add_trace(self):
        self.obj.add_trace()


class ClusterInteractionIndicatorProbability(pl.variables.Beta):
    '''This is the posterior for the probability of a cluster being on
    '''
    def __init__(self, prior, **kwargs):
        '''Parameters

        prior (pl.variables.Beta)
            - prior probability
        **kwargs
            - Other options like graph, value
        '''
        kwargs['name'] = STRNAMES.INDICATOR_PROB
        pl.variables.Beta.__init__(self, a=prior.a.value, b=prior.b.value,
            dtype=float, **kwargs)
        self.add_prior(prior)

    def initialize(self, value_option, hyperparam_option, a=None, b=None, value=None,
        N='auto', delay=0):
        '''Initialize the hyperparameters of the beta prior

        Parameters
        ----------
        value_option : str
            - Option to initialize the value by
            - Options
                - 'manual'
                    - Set the values manually, `value` must be specified
                - 'auto'
                    - Set to the mean of the prior
        hyperparam_option : str
            - If it is a string, then set it by the designated option
            - Options
                - 'manual'
                    - Set the value manually. `a` and `b` must also be specified
                - 'weak-agnostic'
                    - a=b=0.5
                - 'strong-dense'
                    - a = N(N-1), N are the expected number of clusters
                    - b = 0.5
                - 'strong-sparse'
                    - a = 0.5
                    - b = N(N-1), N are the expected number of clusters
                - 'very-strong-sparse'
                    - a = 0.5
                    - b = n_asvs * (n_asvs-1)
        N : str, int
            This is the number of clusters to set the hyperparam options to 
            (if they are dependent on the number of cluster). If 'auto', set to the expected number
            of clusters from a dirichlet process. Else use this number (must be an int).
        a, b : int, float
            - User specified values
            - Only necessary if `hyperparam_option` == 'manual'
        '''
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        self.delay = delay

        if hyperparam_option == 'manual':
            if pl.isnumeric(a) and pl.isnumeric(b):
                self.prior.a.override_value(a)
                self.prior.b.override_value(b)
            else:
                raise ValueError('a ({}) and b ({}) must be numerics (float, int)'.format(
                    a.__class__, b.__class__))
        elif hyperparam_option in ['weak-agnostic', 'auto']:
            self.prior.a.override_value(0.5)
            self.prior.b.override_value(0.5)
        elif hyperparam_option == 'strong-dense':
            if pl.isstr(N):
                if N == 'auto':
                    N = expected_n_clusters(G=self.G)
                else:
                    raise ValueError('`N` ({}) nto recognized'.format(N))
            elif pl.isint(N):
                if N < 0:
                    raise ValueError('`N` ({}) must be positive'.format(N))
            else:
                raise TypeError('`N` ({}) type not recognized'.format(type(N)))
            self.prior.a.override_value(N * (N - 1))
            self.prior.b.override_value(0.5)
        elif hyperparam_option == 'strong-sparse':
            if pl.isstr(N):
                if N == 'auto':
                    N = expected_n_clusters(G=self.G)
                else:
                    raise ValueError('`N` ({}) nto recognized'.format(N))
            elif pl.isint(N):
                if N < 0:
                    raise ValueError('`N` ({}) must be positive'.format(N))
            else:
                raise TypeError('`N` ({}) type not recognized'.format(type(N)))
            self.prior.a.override_value(0.5)
            self.prior.b.override_value((N * (N - 1)))
        elif hyperparam_option == 'very-strong-sparse':
            N = self.G.data.n_asvs
            self.prior.a.override_value(0.5)
            self.prior.b.override_value((N * (N - 1)))
        else:
            raise ValueError('option `{}` not recognized'.format(hyperparam_option))

        if value_option == 'manual':
            if pl.isnumeric(value):
                self.value = value
            else:
                raise ValueError('`value` ({}) must be a numeric (float,int)'.format(
                    value.__class__))
        elif value_option == 'auto':
            self.value = self.prior.mean()/100000
        else:
            raise ValueError('value option "{}" not recognized for indicator prob'.format(
                value_option))

        self.a.value = self.prior.a.value
        self.b.value = self.prior.b.value
        logging.info('Indicator Probability initialization results:\n' \
            '\tprior a: {}\n' \
            '\tprior b: {}\n' \
            '\tvalue: {}'.format(
                self.prior.a.value, self.prior.b.value, self.value))

    def update(self):
        '''Sample the posterior given the data
        '''
        if self.sample_iter < self.delay:
            return
        self.a.value = self.prior.a.value + \
            self.G[REPRNAMES.CLUSTER_INTERACTION_INDICATOR].num_pos_indicators
        self.b.value = self.prior.b.value + \
            self.G[REPRNAMES.CLUSTER_INTERACTION_INDICATOR].num_neg_indicators
        self.sample()
        return self.value


class ClusterInteractionIndicators(pl.variables.Variable):
    '''This is the posterior of the Indicator variables on the interactions
    between clusters. These clusters are not fixed.
    If `value` is not `None`, then we set that to be the initial indicators
    of the cluster interactions
    '''
    def __init__(self, prior, mp=None, relative=True, **kwargs):
        '''Parameters

        prior : pl.variables.Beta
            This is the prior of the variable
        mp : str, None
            If `None`, then there is no multiprocessing.
            If it is a str, then there are two options:
                'debug': pool is done sequentially and not sent to processors
                'full': pool is done at different processors
        relative : bool
            Whether you update using the relative marginal likelihood or not.
        '''
        if not pl.isbool(relative):
            raise TypeError('`relative` ({}) must be a bool'.format(type(relative)))
        if relative:
            if mp is not None:
                raise ValueError('Multiprocessing is slower for rel. Turn mp off')
            self.update = self.update_relative
        else:
            self.update = self.update_direct

        if mp is not None:
            if not pl.isstr(mp):
                raise TypeError('`mp` ({}) must be a str'.format(type(mp)))
            if mp not in ['full', 'debug']:
                raise ValueError('`mp` ({}) not recognized'.format(mp))

        kwargs['name'] = STRNAMES.CLUSTER_INTERACTION_INDICATOR
        pl.variables.Variable.__init__(self, dtype=bool, **kwargs)
        self.n_asvs = len(self.G.data.asvs)
        self.set_value_shape(shape=(self.n_asvs, self.n_asvs))
        self.add_prior(prior)
        self.clustering = self.G[STRNAMES.CLUSTERING_OBJ]
        self.mp = mp
        self.relative = relative

        # parameters used during update
        self.X = None
        self.y = None
        self.process_prec_matrix = None
        self._strr = 'None'

    def initialize(self, delay=0, run_every_n_iterations=1):
        '''Do nothing, the indicators are set in `ClusterInteractionValue`.

        Parameters
        ----------
        delay : int
            How many iterations to delay starting to update the values
        run_every_n_iterations : int
            Which iteration to run on
        '''
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        if not pl.isint(run_every_n_iterations):
            raise TypeError('`run_every_n_iterations` ({}) must be an int'.format(
                type(run_every_n_iterations)))
        if run_every_n_iterations <= 0:
            raise ValueError('`run_every_n_iterations` ({}) must be > 0'.format(
                run_every_n_iterations))

        self.delay = delay
        self.run_every_n_iterations = run_every_n_iterations
        self._there_are_perturbations = self.G.perturbations is not None
        self.update_cnt_indicators()
        self.interactions = self.G[REPRNAMES.INTERACTIONS_OBJ]
        self.n_asvs = len(self.G.data.asvs)

        # These are for the function `self._make_idx_for_clusters`
        self.ndts_bias = []
        self.n_replicates = self.G.data.n_replicates
        self.n_dts_for_replicate = self.G.data.n_dts_for_replicate
        self.total_dts = np.sum(self.n_dts_for_replicate)
        self.replicate_bias = np.zeros(self.n_replicates, dtype=int)
        for ridx in range(1, self.n_replicates):
            self.replicate_bias[ridx] = self.replicate_bias[ridx-1] + \
                self.n_asvs * self.n_dts_for_replicate[ridx - 1]
        for ridx in range(self.G.data.n_replicates):
            self.ndts_bias.append(
                np.arange(0, self.G.data.n_dts_for_replicate[ridx] * self.n_asvs, self.n_asvs))

        # Makes a dictionary that maps the asv index to the rows that it the ASV in
        self.oidx2rows = {}
        for oidx in range(self.n_asvs):
            idxs = np.zeros(self.total_dts, dtype=int)
            i = 0
            for ridx in range(self.n_replicates):
                temp = np.arange(0, self.n_dts_for_replicate[ridx] * self.n_asvs, self.n_asvs)
                temp = temp + oidx
                temp = temp + self.replicate_bias[ridx]
                l = len(temp)
                idxs[i:i+l] = temp
                i += l
            self.oidx2rows[oidx] = idxs

    def add_trace(self):
        self.value = self.G[REPRNAMES.INTERACTIONS_OBJ].get_datalevel_indicator_matrix()
        pl.variables.Variable.add_trace(self)

    def update_cnt_indicators(self):
        self.num_pos_indicators = self.G[REPRNAMES.INTERACTIONS_OBJ].num_pos_indicators()
        self.num_neg_indicators = self.G[REPRNAMES.INTERACTIONS_OBJ].num_neg_indicators()

    def __str__(self):
        return self._strr

    # @profile
    def update_direct(self):
        '''Permute the order that the indices that are updated.

        Build the full master interaction matrix that we can then slice
        '''
        start = time.time()
        if self.sample_iter < self.delay:
            # for interaction in self.interactions:
            #     interaction.indicator=False
            self._strr = '{}\ntotal time: {}'.format(
                self.interactions.get_indicators(), time.time()-start)
            return
        if self.sample_iter % self.run_every_n_iterations != 0:
            return

        # keys = npr.permutation(self.interactions.key_pairs())
        idxs = npr.permutation(self.interactions.size)
        for idx in idxs:
            # print('indicator {}/{}'.format(iii, len(keys)))
            self.update_single_idx_slow(idx=idx)

        self.update_cnt_indicators()
        # Since slicing is literally so slow, it is faster to build than just slicing M
        self.G.data.design_matrices[REPRNAMES.CLUSTER_INTERACTION_VALUE].M.build(
            build=True, build_for_neg_ind=False)
        iii = self.interactions.get_indicators()
        n_on = np.sum(iii)
        self._strr = '{}\ntotal time: {}, n_interactions: {}/{}, {:.2f}'.format(
            iii, time.time()-start, n_on, len(iii), n_on/len(iii))

    def update_single_idx_slow(self, idx):
        '''Update the likelihood for interaction `idx`

        Parameters
        ----------
        idx : int
            This is the index of the interaction we are updating
        '''
        prior_ll_on = np.log(self.G[REPRNAMES.INDICATOR_PROB].value)
        prior_ll_off = np.log(1 - self.G[REPRNAMES.INDICATOR_PROB].value)

        d_on = self.calculate_marginal_loglikelihood(idx=idx, val=True)
        d_off = self.calculate_marginal_loglikelihood(idx=idx, val=False)

        ll_on = d_on['ret'] + prior_ll_on
        ll_off = d_off['ret'] + prior_ll_off

        # print('slow\n\ttotal: {}\n\tbeta_logdet_diff: {}\n\t' \
        #     'priorvar_logdet_diff: {}\n\tbEb_diff: {}\n\t' \
        #     'bEbprior_diff: {}\n\tn_on_when_off: {}'.format(
        #         ll_on - ll_off,
        #         d_on['beta_logdet'] - d_off['beta_logdet'],
        #         d_on['priorvar_logdet'] - d_off['priorvar_logdet'],
        #         d_on['bEb'] - d_off['bEb'],
        #         d_on['bEbprior'] - d_off['bEbprior'],
        #         self.interactions.num_pos_indicators()))

        dd = [ll_off, ll_on]

        res = bool(sample_categorical_log(dd))
        self.interactions.iloc(idx).indicator = res
        self.update_cnt_indicators()

    # @profile
    def _make_idx_vector_for_clusters(self):
        '''Creates a dictionary that maps the cluster id to the
        rows that correspond to each ASV in the cluster.

        We cannot cast this with numba because it does not support Fortran style
        raveling :(.

        Returns
        -------
        dict: int -> np.ndarray
            Maps the cluster ID to the row indices corresponding to it
        '''
        clusters = [np.asarray(oidxs, dtype=int).reshape(-1,1) \
            for oidxs in self.clustering.tolistoflists()]

        d = {}
        cids = self.clustering.order

        for cidx,cid in enumerate(cids):
            a = np.zeros(len(clusters[cidx]) * self.total_dts, dtype=int)
            i = 0
            for ridx in range(self.n_replicates):
                idxs = np.zeros(
                    (len(clusters[cidx]),
                    self.n_dts_for_replicate[ridx]), int)
                idxs = idxs + clusters[cidx]
                idxs = idxs + self.ndts_bias[ridx]
                idxs = idxs + self.replicate_bias[ridx]
                idxs = idxs.ravel('F')
                l = len(idxs)
                a[i:i+l] = idxs
                i += l

            d[cid] = a
        
        if self.G.data.zero_inflation_transition_policy is not None:
            # We need to convert the indices that are meant from no zero inflation to 
            # ones that take into account zero inflation - use the array from 
            # `data.Data._setrows_to_include_zero_inflation`. If the index should be
            # included, then we subtract the number of indexes that are previously off
            # before that index. If it should not be included then we exclude it
            prevoff_arr = self.G.data.off_previously_arr_zero_inflation
            rows_to_include = self.G.data.rows_to_include_zero_inflation
            for cid in d:
                arr = d[cid]
                new_arr = np.zeros(len(arr), dtype=int)
                n = 0
                for idx in arr:
                    if rows_to_include[idx]:
                        new_arr[n] = idx - prevoff_arr[idx]
                        n += 1

                new_arr = new_arr[:n]
                d[cid] = new_arr
        return d

    # @profile
    def make_rel_params(self):
        '''We make the parameters needed to update the relative log-likelihod.
        This function is called once at the beginning of the update.

        Parameters that we create with this function
        --------------------------------------------
        ys : dict (int -> np.ndarray)
            Maps the target cluster id to the observation matrix that it
            corresponds to (only the ASVs in the target cluster). This 
            array already has the growth and self-interactions subtracted
            out:
                $ \frac{log(x_{k+1}) - log(x_{k})}{dt} - a_{1,k} - a_{2,k}x_{k} $
        process_precs : dict (int -> np.ndarray)
            Maps the target cluster id to the vector of the process precision
            that corresponds to the target cluster (only the ASVs in the target
            cluster). This is a 1D array that corresponds to the diagonal of what
            would be the precision matrix.
        interactionXs : dict (int -> np.ndarray)
            Maps the target cluster id to the matrix of the design matrix of the
            interactions. Only includes the rows that correspond to the ASVs in the
            target cluster. It includes every single column as if all of the indicators
            are on. We only index out the columns when we are doing the marginalization.
        prior_prec_interaction : float
            Prior precision of the interaction value. We then use this
            value to make the diagonal of the prior precision.
        prior_var_interaction : float
            Prior variance of the interaction value.
        prior_mean_interaction : float
            Prior mean of the interaction values. We use this value
            to make the prior mean vector during the marginalization.
        n_on_master : int
            How many interactions are on at any one time. We adjust this
            throughout the update depending on what interactions we turn off and
            on.
        prior_ll_on : float
            Prior log likelihood of a positive interaction
        prior_ll_off : float
            Prior log likelihood of the negative interaction
        priorvar_logdet_diff : float
            This is the prior variance log determinant that we add when the indicator
            is positive.

        Parameters created if there are perturbations
        ---------------------------------------------
        perturbationsXs : dict (int -> np.ndarray)
            Maps the target cluster id to the design matrix that corresponds to 
            the on perturbations of the target clusters. This is preindexed in
            both rows and columns
        prior_prec_perturbations : dict (int -> np.ndarray)
            Maps the target cluster id to the diagonal of the prior precision
            of the perturbations
        prior_var_perturbations : dict (int -> np.ndarray)
            Maps the target cluster id to the diagonal of the prior variance 
            of the perturbations
        prior_mean_perturbations : dict (int -> np.ndarray)
            Maps the target cluster id to the vector of the prior mean of the
            perturbations
        '''
        # Get the row indices for each cluster
        row_idxs = self._make_idx_vector_for_clusters()

        # Create ys
        self.ys = {}
        y = self.G.data.construct_lhs(keys=[
            REPRNAMES.SELF_INTERACTION_VALUE, REPRNAMES.GROWTH_VALUE],
            kwargs_dict={REPRNAMES.GROWTH_VALUE:{'with_perturbations': False}})
        for tcid in self.clustering.order:
            self.ys[tcid] = y[row_idxs[tcid], :]

        # Create process_precs
        self.process_precs = {}
        process_prec_diag = self.G[REPRNAMES.PROCESSVAR].prec
        for tcid in self.clustering.order:
            self.process_precs[tcid] = process_prec_diag[row_idxs[tcid]]

        # Make interactionXs
        self.interactionXs = {}
        self.G.data.design_matrices[REPRNAMES.CLUSTER_INTERACTION_VALUE].M.build(
            build=True, build_for_neg_ind=True)
        XM_master = self.G.data.design_matrices[REPRNAMES.CLUSTER_INTERACTION_VALUE].toarray()
        for tcid in self.clustering.order:
            self.interactionXs[tcid] = XM_master[row_idxs[tcid], :]

        # Make prior parameters
        self.prior_var_interaction = self.G[REPRNAMES.PRIOR_VAR_INTERACTIONS].value
        self.prior_prec_interaction = 1/self.prior_var_interaction
        self.prior_mean_interaction = self.G[REPRNAMES.PRIOR_MEAN_INTERACTIONS].value
        self.prior_ll_on = np.log(self.prior.value)
        self.prior_ll_off = np.log(1 - self.prior.value)
        self.n_on_master = self.interactions.num_pos_indicators()

        # Make priorvar_logdet
        self.priorvar_logdet = np.log(self.prior_var_interaction)

        if self._there_are_perturbations:
            XMpert_master = self.G.data.design_matrices[REPRNAMES.PERT_VALUE].toarray()

            # Make perturbationsXs
            self.perturbationsXs = {}
            for tcid in self.clustering.order:
                rows = row_idxs[tcid]
                cols = []
                i = 0
                for perturbation in self.G.perturbations:
                    for cid in perturbation.indicator.value:
                        if perturbation.indicator.value[cid]:
                            if cid == tcid:
                                cols.append(i)
                            i += 1
                cols = np.asarray(cols, dtype=int)

                self.perturbationsXs[tcid] = pl.util.fast_index(M=XMpert_master,
                    rows=rows, cols=cols)

            # Make prior perturbation parameters
            self.prior_mean_perturbations = {}
            self.prior_var_perturbations = {}
            self.prior_prec_perturbations = {}
            for tcid in self.clustering.order:
                mean = []
                var = []
                for perturbation in self.G.perturbations:
                    if perturbation.indicator.value[tcid]:
                        # This is on, get the parameters
                        mean.append(perturbation.magnitude.prior.mean.value)
                        var.append(perturbation.magnitude.prior.var.value)
                self.prior_mean_perturbations[tcid] = np.asarray(mean)
                self.prior_var_perturbations[tcid] = np.asarray(var)
                self.prior_prec_perturbations[tcid] = 1/self.prior_var_perturbations[tcid]

            # Make priorvar_det_perturbations
            self.priorvar_det_perturbations = 0
            for perturbation in self.G.perturbations:
                self.priorvar_det_perturbations += \
                    perturbation.indicator.num_on_clusters() * \
                    perturbation.magnitude.prior.var.value

    # @profile
    def update_relative(self):
        '''Update the indicators variables by calculating the relative loglikelihoods
        of it being on as supposed to off. Because this is a relative loglikelihood,
        we only need to take into account the following parameters of the model:
            - Only the ASVs in the target cluster of the interaction
            - Only the positively indicated interactions going into the
              target cluster.

        This is 1000's of times faster than `update` because we are operating on matrices
        that are MUCH smaller than in a full system. These matrices are also considered dense
        so we do all of our computations without sparse matrices.

        We permute the order that the indices are updated for more robust mixing.
        '''
        start = time.time()
        if self.sample_iter < self.delay:
            self._strr = '{}\ntotal time: {}'.format(
                self.interactions.get_indicators(), time.time()-start)
            return
        if self.sample_iter % self.run_every_n_iterations != 0:
            return

        idxs = npr.permutation(self.interactions.size)

        self.make_rel_params()
        for idx in idxs:
            self.update_single_idx_fast(idx=idx)

        self.update_cnt_indicators()
        # Since slicing is literally so slow, it is faster to build than just slicing M
        self.G.data.design_matrices[REPRNAMES.CLUSTER_INTERACTION_VALUE].M.build(
            build=True, build_for_neg_ind=False)
        iii = self.interactions.get_indicators()
        n_on = np.sum(iii)
        self._strr = '{}\ntotal time: {}, n_interactions: {}/{}, {:.2f}'.format(
            iii, time.time()-start, n_on, len(iii), n_on/len(iii))

    def calculate_marginal_loglikelihood(self, idx, val):
        '''Calculate the likelihood of interaction `idx` with the value `val`
        '''
        # Build and initialize
        self.interactions.iloc(idx).indicator = val
        self.update_cnt_indicators()
        self.G.data.design_matrices[REPRNAMES.CLUSTER_INTERACTION_VALUE].M.build()

        lhs = [REPRNAMES.GROWTH_VALUE, REPRNAMES.SELF_INTERACTION_VALUE]
        if self._there_are_perturbations:
            rhs = [REPRNAMES.PERT_VALUE, REPRNAMES.CLUSTER_INTERACTION_VALUE]
        else:
            rhs = [REPRNAMES.CLUSTER_INTERACTION_VALUE]

        y = self.G.data.construct_lhs(lhs, 
            kwargs_dict={REPRNAMES.GROWTH_VALUE:{'with_perturbations': False}})
        X = self.G.data.construct_rhs(rhs, toarray=True)

        if X.shape[1] == 0:
            return {
            'ret': 0,
            'beta_logdet': 0,
            'priorvar_logdet': 0,
            'bEb': 0,
            'bEbprior': 0}

        process_prec = self.G[REPRNAMES.PROCESSVAR].build_matrix(cov=False, sparse=False)
        prior_prec = build_prior_covariance(G=self.G, cov=False, order=rhs, sparse=False)
        prior_var = build_prior_covariance(G=self.G, cov=True, order=rhs, sparse=False)
        prior_mean = build_prior_mean(G=self.G, order=rhs, shape=(-1,1))


        # Calculate the posterior
        beta_prec = X.T @ process_prec @ X + prior_prec
        beta_cov = pinv(beta_prec, self)
        beta_mean = beta_cov @ ( X.T @ process_prec @ y + prior_prec @ prior_mean )
        
        beta_mean = np.asarray(beta_mean).reshape(-1,1)

        # Perform the marginalization
        try:
            beta_logdet = log_det(beta_cov, self)
        except:
            logging.critical('Crashed in log_det')
            logging.critical('beta_cov:\n{}'.format(beta_cov))
            logging.critical('prior_prec\n{}'.format(prior_prec))
            raise
        priorvar_logdet = log_det(prior_var, self)
        ll2 = 0.5 * (beta_logdet - priorvar_logdet)

        a = np.asarray(prior_mean.T @ prior_prec @ prior_mean)[0,0]
        b = np.asarray(beta_mean.T @ beta_prec @ beta_mean)[0,0]
        ll3 = -0.5 * (a  - b)

        return {
            'ret': ll2+ll3,
            'beta_logdet': beta_logdet,
            'priorvar_logdet': priorvar_logdet,
            'bEb': b,
            'bEbprior': a}

    def update_single_idx_fast(self, idx):
        '''Calculate the relative log likelihood of changing the indicator of the
        interaction at index `idx`.

        This is about 20X faster than `update_single_idx_slow`.

        Parameters
        ----------
        idx : int
            This is the index of the interaction index we are sampling
        '''
        # Get the current interaction by the index
        self.curr_interaction = self.interactions.iloc(idx)
        start_sign = self.curr_interaction.indicator
        
        tcid = self.curr_interaction.target_cid
        self.curr_interaction.indicator = False
        self.col_idxs = np.asarray(
            self.interactions.get_arg_indicators(target_cid=tcid),
            dtype=int)

        if not start_sign:
            self.n_on_master += 1
            self.num_pos_indicators += 1
            self.num_neg_indicators -= 1

        d_on = self.calculate_relative_marginal_loglikelihood(idx=idx, val=True)

        self.n_on_master -= 1
        self.num_pos_indicators -= 1
        self.num_neg_indicators += 1
        d_off = self.calculate_relative_marginal_loglikelihood(idx=idx, val=False)

        ll_on = d_on + self.prior_ll_on
        ll_off = d_off + self.prior_ll_off

        dd = [ll_off, ll_on]

        # print('\nindicator', idx)
        # print('fast\n\ttotal: {}\n\tbeta_logdet_diff: {}\n\t' \
        #     'priorvar_logdet_diff: {}\n\tbEb_diff: {}\n\t' \
        #     'bEbprior_diff: {}\n\tn_on_when_off: {}'.format(
        #         ll_on - ll_off,
        #         d_on['beta_logdet'] - d_off['beta_logdet'],
        #         d_on['priorvar_logdet'] - d_off['priorvar_logdet'],
        #         d_on['bEb'] - d_off['bEb'],
        #         d_on['bEbprior'] - d_off['bEbprior'],
        #         self.n_on_master))
        # print('log(prior_var)', np.log(self.prior_var_interaction))
        # self.update_single_idx_slow(idx)

        res = bool(sample_categorical_log(dd))
        if res:
            self.n_on_master += 1
            self.num_pos_indicators += 1
            self.num_neg_indicators -= 1
            self.curr_interaction.indicator = True

    def calculate_relative_marginal_loglikelihood(self, idx, val):
        '''Calculate the relative marginal log likelihood for the interaction index
        `idx` with the indicator `val`

        Parameters
        ----------
        idx : int
            This is the index of the interaction
        val : bool
            This is the value to calculate it as
        '''
        tcid = self.curr_interaction.target_cid

        y = self.ys[tcid]
        process_prec = self.process_precs[tcid]

        # Make X, prior mean, and prior_var
        if val:
            cols = np.append(self.col_idxs, idx)
        else:
            cols = self.col_idxs
        
        X = self.interactionXs[tcid][:, cols]
        prior_mean = np.full(len(cols), self.prior_mean_interaction)
        prior_prec_diag = np.full(len(cols), self.prior_prec_interaction)

        if self._there_are_perturbations:
            Xpert = self.perturbationsXs[tcid]
            X = np.hstack((X, Xpert))

            prior_mean = np.append(
                prior_mean,
                self.prior_mean_perturbations[tcid])
            
            prior_prec_diag = np.append(
                prior_prec_diag,
                self.prior_prec_perturbations[tcid])

        if X.shape[1] == 0:
            # return {
            #     'ret': 0,
            #     'beta_logdet': 0,
            #     'priorvar_logdet': 0,
            #     'bEb': 0,
            #     'bEbprior': 0}
            return 0

        prior_prec = np.diag(prior_prec_diag)
        pm = (prior_prec_diag * prior_mean).reshape(-1,1)

        # Do the marginalization
        a = X.T * process_prec

        beta_prec = (a @ X) + prior_prec
        beta_cov = pinv(beta_prec, self)
        beta_mean = beta_cov @ ((a @ y) + pm )

        bEb = (beta_mean.T @ beta_prec @ beta_mean)[0,0]
        try:
            beta_logdet = log_det(beta_cov, self)
        except:
            logging.critical('Crashed in log_det')
            logging.critical('beta_cov:\n{}'.format(beta_cov))
            logging.critical('prior_prec\n{}'.format(prior_prec))

            logging.critical('here')
            print('y')
            print(y.shape)
            print(y)
            print('process_prec')
            print(process_prec.shape)
            print('X')
            print(X.shape)
            print(X)
            print('priors')
            print(prior_mean)
            print(prior_prec_diag)
            print('self-interactions')
            X = pl.toarray(self.G.data.design_matrices[REPRNAMES.SELF_INTERACTION_VALUE].matrix)
            print(X.shape)
            print(np.any(np.isnan(X)))
            print('growth')
            X = pl.toarray(self.G.data.design_matrices[REPRNAMES.GROWTH_VALUE].matrix_without_perturbations)
            print(X.shape)
            print(np.any(np.isnan(X)))
            print('orig y')
            y = self.G.data.lhs.vector
            print(y.shape)
            print(np.any(np.isnan(y)))
            print('cluster-interactions')
            X = pl.toarray(self.G.data.design_matrices[REPRNAMES.CLUSTER_INTERACTION_VALUE].matrix)
            print(X.shape)
            print(np.any(np.isnan(X)))

            n_on = 0
            for row in range(X.shape[0]):
                n_on += np.any(np.isnan(X[row]))

            print('nans on {}/{} rows'.format(n_on, X.shape[0]))
                    

            raise
        
        if val:
            bEbprior = (self.prior_mean_interaction**2)/self.prior_var_interaction
            priorvar_logdet = self.priorvar_logdet
        else:
            bEbprior = 0
            priorvar_logdet = 0

        ll2 = 0.5 * (beta_logdet - priorvar_logdet)
        ll3 = 0.5 * (bEb - bEbprior)

        # return {
        #     'ret': ll2+ll3,
        #     'beta_logdet': beta_logdet,
        #     'priorvar_logdet': priorvar_logdet,
        #     'bEb': bEb,
        #     'bEbprior': bEbprior}
        return ll2 + ll3
            
    def kill(self):
        pass


# Logistic Growth
# ---------------
class PriorVarMH(pl.variables.SICS):
    '''This is the posterior for the prior variance of either the growth
    or self-interaction parameter. We update with a MH update since this 
    prior is not conjugate.

    Parameters
    ----------
    prior : pl.variables.SICS
        This is the prior of this distribution - which is a Squared
        Inverse Chi Squared (SICS) distribution
    child_name : str
        This is the name of the variable that this is a prior variance
        for. This is either the name of the growth parameter or the 
        self-interactions parameter
    kwargs : dict
        These are the other parameters for the initialization.
    '''

    def __init__(self, prior, child_name, **kwargs):
        if child_name == STRNAMES.GROWTH_VALUE:
            kwargs['name'] = STRNAMES.PRIOR_VAR_GROWTH
        elif child_name == STRNAMES.SELF_INTERACTION_VALUE:
            kwargs['name'] = STRNAMES.PRIOR_VAR_SELF_INTERACTIONS
        else:
            raise ValueError('`child_name` ({}) not recognized'.format(child_name))
        pl.variables.SICS.__init__(self, dtype=float, **kwargs)
        self.child_name = child_name
        self.add_prior(prior)
        self.proposal = pl.variables.SICS(dof=None, scale=None, value=None)

    def __str__(self):
        # If this fails, it is because we are dividing by 0 sampler_iter
        # If which case we just return the value 
        try:
            s = 'Value: {}, Acceptance rate: {}'.format(
                self.value, np.mean(self.acceptances[
                    np.max([self.sample_iter-50, 0]):self.sample_iter]))
        except:
            s = str(self.value)
        return s

    def initialize(self, value_option, dof_option, scale_option, 
        proposal_option, target_acceptance_rate, tune, end_tune,
        value=None, dof=None, scale=None, proposal_dof=None, delay=0):
        '''Initialize the parameters of the distribution and the 
        proposal distribution

        Parameters
        ----------
        value_option : str
            Different ways to initialize the values
            Options
                'manual'
                    Set the value manually, `value` must also be specified
                'unregularized'
                    Do unregularized regression and the value is set to the
                    variance of the growth values
                'prior-mean', 'auto'
                    Set the value to the prior of the mean
        scale_option : str
            Different ways to initialize the scale of the prior
            Options
                'manual'
                    Set the value manually, `scale` must also be specified
                'auto', 'inflated-median'
                    We set the scale such that the mean of the prior is
                    equal to the median growth values calculated
                    with linear regression squared and inflated by 100.
        dof_option : str
            How informative the prior should be (setting the dof)
                'diffuse': set to the mimumum value (2)
                'weak': set so that 10% of the posterior comes from the prior
                'strong': set so that 50% of the posterior comes from the prior
                'manual': set to the value provided in the parameter `shape`
                'auto': Set to 'weak'
        proposal_option : str
            How to set the initial dof of the proposal - this will get adjusted with
            tuning
                'tight', 'auto'
                    Set the dof to be 15, relatively strong initially
                'diffuse'
                    Set the dof to be 2.5, relatively diffuse initially
                'manual'
                    Set the dof with the parameter `proposal_dof'
        target_acceptance_rate : float, str
            This is the target_acceptance rate. Options:
                'auto', 'optimal'
                    Set to 0.44
                float
                    This is the value you want
        tune : str, int
            This is how often you want to update the proposal dof
                int
                'auto'
                    Set to every 50 iterations
        end_tune : str, int
            This is when to stop the tuning
                'half-burnin', 'auto'
                    Half of burnin, rounded down
                int
        delay : int
            How many iterations to delay updating the value of the variance
        '''
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        self.delay = delay

        # Set the proposal dof
        if not pl.isstr(proposal_option):
            raise TypeError('`proposal_option` ({}) must be a str'.format(
                type(proposal_option)))
        elif proposal_option == 'manual':
            if not pl.isnumeric(proposal_dof):
                raise TypeError('`proposal_dof` ({}) must be a numeric'.format(
                    type(proposal_dof)))
            if proposal_dof < 2:
                raise ValueError('`proposal_dof` ({}) not proper'.format(proposal_dof))
        elif proposal_option in ['tight', 'auto']:
            proposal_dof = 15
        elif proposal_option == 'diffuse':
            proposal_dof = 2.5
        else:
            raise ValueError('`proposal_option` ({}) not recognized'.format(
                proposal_option))
        self.proposal.dof.value = proposal_dof

        # Set the propsal parameters
        if pl.isstr(target_acceptance_rate):
            if target_acceptance_rate in ['optimal', 'auto']:
                target_acceptance_rate = 0.44
            else:
                raise ValueError('`target_acceptance_rate` ({}) not recognized'.format(
                    target_acceptance_rate))
        elif pl.isfloat(target_acceptance_rate):
            if target_acceptance_rate < 0 or target_acceptance_rate > 1:
                raise ValueError('`target_acceptance_rate` ({}) out of range'.format(
                    target_acceptance_rate))
        else:
            raise TypeError('`target_acceptance_rate` ({}) type not recognized'.format(
                type(target_acceptance_rate)))
        self.target_acceptance_rate = target_acceptance_rate

        if pl.isstr(tune):
            if tune in ['auto']:
                tune = 50
            else:
                raise ValueError('`tune` ({}) not recognized'.format(tune))
        elif pl.isint(tune):
            if tune < 0:
                raise ValueError('`tune` ({}) must be > 0'.format(
                    tune))
        else:
            raise TypeError('`tune` ({}) type not recognized'.format(type(tune)))
        self.tune = tune

        if pl.isstr(end_tune):
            if end_tune in ['auto', 'half-burnin']:
                end_tune = int(self.G.inference.burnin/2)
            else:
                raise ValueError('`tune` ({}) not recognized'.format(end_tune))
        elif pl.isint(end_tune):
            if end_tune < 0 or end_tune > self.G.inference.burnin:
                raise ValueError('`end_tune` ({}) out of range (0, {})'.format(
                    end_tune, self.G.inference.burnin))
        else:
            raise TypeError('`end_tune` ({}) type not recognized'.format(type(end_tune)))
        self.end_tune = end_tune

        # Set the prior dof
        if not pl.isstr(dof_option):
            raise TypeError('`dof_option` ({}) must be a str'.format(type(dof_option)))
        if dof_option == 'manual':
            if not pl.isnumeric(dof):
                raise TypeError('`dof` ({}) must be a numeric'.format(type(dof)))
            if dof < 2:
                raise ValueError('`dof` ({}) must be >= 2'.format(dof))
        elif dof_option == 'diffuse':
            dof = 2.5
        elif dof_option in ['weak', 'auto']:
            dof = len(self.G.data.asvs)/9
        elif dof_option == 'strong':
            dof = len(self.G.data.asvs)/2
        else:
            raise ValueError('`dof_option` ({}) not recognized'.format(dof_option))
        if dof < 2:
            raise ValueError('`dof` ({}) must be strictly larger than 2 to be a proper' \
                ' prior'.format(dof))
        self.prior.dof.override_value(dof)

        # Set the prior scale
        if not pl.isstr(scale_option):
            raise TypeError('`scale_option` ({}) must be a str'.format(type(scale_option)))
        if scale_option == 'manual':
            if not pl.isnumeric(scale):
                raise TypeError('`scale` ({}) must be a numeric'.format(type(scale)))
            if scale <= 0:
                raise ValueError('`scale` ({}) must be positive'.format(scale))
        elif scale_option in ['auto', 'inflated-median']:
            # Perform linear regression
            rhs = [REPRNAMES.GROWTH_VALUE, REPRNAMES.SELF_INTERACTION_VALUE]
            X = self.G.data.construct_rhs(keys=rhs,
                kwargs_dict={REPRNAMES.GROWTH_VALUE:{'with_perturbations':False}},
                index_out_perturbations=True)
            y = self.G.data.construct_lhs(index_out_perturbations=True)

            prec = X.T @ X
            cov = pinv(prec, self)
            mean = cov @ X.T @ y
            if self.child_name == STRNAMES.GROWTH_VALUE:
                mean = 1e4*(np.median(mean[:self.G.data.n_asvs]) ** 2)
            else:
                mean = 1e4*(np.median(mean[self.G.data.n_asvs:]) ** 2)

            # Calculate the scale
            scale = mean * (self.prior.dof.value - 2) / self.prior.dof.value
        else:
            raise ValueError('`scale_option` ({}) not recognized'.format(scale_option))
        self.prior.scale.override_value(scale)

        # Set the initial value of the prior
        if not pl.isstr(value_option):
            raise TypeError('`value_option` ({}) must be a str'.format(type(value_option)))
        if value_option == 'manual':
            if not pl.isnumeric(value):
                raise ValueError('If `value_option` == "manual", value ({}) ' \
                    'must be a numeric (float, int)'.format(value.__class__))
        elif value_option in ['inflated-median']:
            # No interactions
            rhs = [
                REPRNAMES.GROWTH_VALUE,
                REPRNAMES.SELF_INTERACTION_VALUE]
            X = self.G.data.construct_rhs(keys=rhs,
                kwargs_dict={REPRNAMES.GROWTH_VALUE:{'with_perturbations':False}},
                index_out_perturbations=True)
            y = self.G.data.construct_lhs(index_out_perturbations=True)

            prec = X.T @ X
            cov = pinv(prec, self)
            mean = cov @ X.T @ y
            if self.child_name == STRNAMES.GROWTH_VALUE:
                value = 1e4*(np.median(mean[:self.G.data.n_asvs]) ** 2)
            else:
                value = 1e4*(np.median(mean[self.G.data.n_asvs:]) ** 2)
        elif value_option in ['prior-mean', 'auto']:
            value = self.prior.mean()
        else:
            raise ValueError('`value_option` "{}" not recognized'.format(value_option))
        self.value = value

    def update_dof(self):
        '''Updat the `dof` parameter so that we adjust the acceptance
        rate to `target_acceptance_rate`
        '''
        if self.sample_iter == 0:
            self.temp_acceptances = 0
            self.acceptances = np.zeros(self.G.inference.n_samples, dtype=bool)
        
        elif self.sample_iter > self.end_tune:
            # Don't do any more updates
            return
        
        elif self.sample_iter % self.tune == 0:
            # Update dof
            acceptance_rate = self.temp_acceptances / self.tune
            if acceptance_rate > self.target_acceptance_rate:
                self.proposal.dof.value = self.proposal.dof.value * 1.5
            else:
                self.proposal.dof.value = self.proposal.dof.value / 1.5
            self.temp_acceptances = 0

    def update(self):
        '''First we check if we need to tune the dof, which we do during
        the first half of burnin. We calculate the likelihoods in logspace
        '''
        if self.sample_iter < self.delay:
            return
        self.update_dof()

        # Get necessary data of the respective parameter
        var = self.G[self.child_name]
        x = var.value.ravel()
        mu = var.prior.mean.value
        low = var.low
        high = var.high

        # propose a new value
        prev_value = self.value
        prev_value_std = math.sqrt(prev_value)
        self.proposal.scale.value = self.value
        new_value = self.proposal.sample() # Sample a new value
        new_value_std = math.sqrt(new_value)

        # Calculate the target distribution ll
        prev_target_ll = 0
        for i in range(len(x)):
            prev_target_ll += pl.random.truncnormal.logpdf(
                value=x[i], mean=mu, std=prev_value_std,
                low=low, high=high)
        new_target_ll = 0
        for i in range(len(x)):
            new_target_ll += pl.random.truncnormal.logpdf(
                value=x[i], mean=mu, std=new_value_std,
                low=low, high=high)

        # Normalize by the ll of the proposal
        prev_prop_ll = self.proposal.logpdf(value=prev_value)
        new_prop_ll = self.proposal.logpdf(value=new_value)

        # Accept or reject
        r = (new_target_ll - prev_prop_ll) - \
            (prev_target_ll - new_prop_ll)
        u = np.log(pl.random.misc.fast_sample_standard_uniform())

        # print('\n\n\n{} prior_var\n----------'.format(self.child_name))
        # print('prev_value', prev_value)
        # print('prev_target_ll', prev_target_ll)
        # print('prev_prop_ll', prev_prop_ll)
        # print('new value', new_value)
        # print('new_target_ll', new_target_ll)
        # print('new_prop_ll', new_prop_ll)
        # print('mu', mu)
        # print('prev_value_std', prev_value_std)
        # print('new_value_std', new_value_std)
        # print('\nr', r, u)

        if r >= u:
            self.acceptances[self.sample_iter] = True
            self.value = new_value
            self.temp_acceptances += 1
        else:
            self.value = prev_value

    def visualize(self, path, f, section='posterior'):
        '''Render the traces in the folder `basepath` and write the 
        learned values to the file `f`.

        Parameters
        ----------
        path : str
            This is the path to write the files to
        f : _io.TextIOWrapper
            File that we are writing the values to
        section : str
            Section of the trace to compute on. Options:
                'posterior' : posterior samples
                'burnin' : burn-in samples
                'entire' : both burn-in and posterior samples

        Returns
        -------
        _io.TextIOWrapper
        '''
        f.write('\n\n###################################\n')
        f.write(self.name)
        f.write('\n###################################\n')
        if not self.G.inference.is_being_traced(self):
            f.write('`{}` not learned\n\tValue: {}\n'.format(self.name, self.value))
            return f
        summ = pl.summary(self, section=section)
        for k,v in summ.items():
            f.write('\t{}: {}\n'.format(k,v))

        # Plot the traces
        ax1, ax2 = visualization.render_trace(var=self, plt_type='both', section=section,
            include_burnin=True, log_scale=True, rasterized=True)

        # Plot the prior over the posterior
        l,h = ax1.get_xlim()
        xs = np.arange(l,h,step=(h-l)/100) 
        ys = []
        for x in xs:
            ys.append(pl.random.sics.pdf(value=x, 
                dof=self.prior.dof.value,
                scale=self.prior.scale.value))
        ax1.plot(xs, ys, label='prior', alpha=0.5, color='red', rasterized=True)
        ax1.legend()

        # Plot the acceptance rate over the trace
        ax3 = ax2.twinx()
        ax3 = visualization.render_acceptance_rate_trace(var=self, ax=ax3, 
            label='Acceptance Rate', color='red', scatter=False, rasterized=True)
        ax3.legend()
        fig = plt.gcf()
        fig.tight_layout()
        fig.suptitle(self.name)
        plt.savefig(path)
        plt.close()

        return f


class PriorMeanMH(pl.variables.TruncatedNormal):
    '''This implements the posterior for the prior mean of the either
    the growths or the self-interactions

    Parameters
    ----------
    '''
    def __init__(self, prior, child_name, **kwargs):
        if child_name == STRNAMES.GROWTH_VALUE:
            kwargs['name'] = STRNAMES.PRIOR_MEAN_GROWTH
        elif child_name == STRNAMES.SELF_INTERACTION_VALUE:
            kwargs['name'] = STRNAMES.PRIOR_MEAN_SELF_INTERACTIONS
        else:
            raise ValueError('`child_name` ({}) not recognized'.format(child_name))
        pl.variables.TruncatedNormal.__init__(self, mean=None, var=None, dtype=float, **kwargs)
        self.child_name = child_name
        self.add_prior(prior)
        self.proposal = pl.variables.TruncatedNormal(mean=None, var=None, value=None)

    def __str__(self):
        # If this fails, it is because we are dividing by 0 sampler_iter
        # If which case we just return the value 
        try:
            s = 'Value: {}, Acceptance rate: {}'.format(
                self.value, np.mean(self.acceptances[
                    np.max([self.sample_iter-50, 0]):self.sample_iter]))
        except:
            s = str(self.value)
        return s

    def initialize(self, value_option, mean_option, var_option,
        truncation_settings, proposal_option, target_acceptance_rate,
        tune, end_tune, value=None, mean=None, var=None, proposal_var=None,
        delay=0):
        '''These are the parameters to initialize the parameters
        of the class. Depending whether it is a self-interaction
        or a growth, it does it differently.

        Parameters
        ----------
        value_option : str
            How to initialize the value. Options:
                'auto', 'prior-mean'
                    Set to the prior mean
                'linear-regression'
                    Set the values from an unregularized linear regression
                'manual'
                    `value` must also be specified
        truncation_settings: str, tuple
            How to set the truncation parameters. The proposal trucation will
            be set the same way.
                tuple - (low,high)
                    These are the truncation parameters
                'auto'
                    If self-interactions, 'negative'. If growths, 'positive'
                'positive'
                    (0, \infty)
                'negative'
                    (-\infty, 0)
                'in-vivo'
                    Not implemented
        mean_option : str
            How to set the mean
                'auto', 'median-linear-regression'
                    Set the mean to the median of the values from an
                    unregularized linear-regression
                'manual'
                    `mean` must also be specified
        var_option : str
            How to set the var
                'auto', 'diffuse-linear-regression'
                    Set the var to 10^4 * median(a_l)
                'manaul'
                    `var` must also be specified.
        proposal_option : str
            How to initialize the proposal variance:
                'auto'
                    mean**2 / 100
                'manual'
                    `proposal_var` must also be supplied
        target_acceptance_rate : str, float
            If float, this is the target acceptance rate
            If str: 
                'optimal', 'auto': 0.44
        tune : str, int
            How often to tune the proposal. If str:
                'auto': 50
        end_tune : str, int
            When to stop tuning the proposal. If str:
                'auto', 'half-burnin': Half of burnin
        '''
        self._there_are_perturbations = self.G.perturbations is not None
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        self.delay = delay

        # Set the propsal parameters
        if pl.isstr(target_acceptance_rate):
            if target_acceptance_rate in ['optimal', 'auto']:
                target_acceptance_rate = 0.44
            else:
                raise ValueError('`target_acceptance_rate` ({}) not recognized'.format(
                    target_acceptance_rate))
        elif pl.isfloat(target_acceptance_rate):
            if target_acceptance_rate < 0 or target_acceptance_rate > 1:
                raise ValueError('`target_acceptance_rate` ({}) out of range'.format(
                    target_acceptance_rate))
        else:
            raise TypeError('`target_acceptance_rate` ({}) type not recognized'.format(
                type(target_acceptance_rate)))
        self.target_acceptance_rate = target_acceptance_rate

        if pl.isstr(tune):
            if tune in ['auto']:
                tune = 50
            else:
                raise ValueError('`tune` ({}) not recognized'.format(tune))
        elif pl.isint(tune):
            if tune < 0:
                raise ValueError('`tune` ({}) must be > 0'.format(
                    tune))
        else:
            raise TypeError('`tune` ({}) type not recognized'.format(type(tune)))
        self.tune = tune

        if pl.isstr(end_tune):
            if end_tune in ['auto', 'half-burnin']:
                end_tune = int(self.G.inference.burnin/2)
            else:
                raise ValueError('`tune` ({}) not recognized'.format(end_tune))
        elif pl.isint(end_tune):
            if end_tune < 0 or end_tune > self.G.inference.burnin:
                raise ValueError('`end_tune` ({}) out of range (0, {})'.format(
                    end_tune, self.G.inference.burnin))
        else:
            raise TypeError('`end_tune` ({}) type not recognized'.format(type(end_tune)))
        self.end_tune = end_tune

        # Set the truncation settings
        if truncation_settings is None:
            truncation_settings = 'positive'
        if pl.isstr(truncation_settings):
            if truncation_settings == 'positive':
                self.low = 0.
                self.high = float('inf')
            # elif truncation_settings == 'negative':
            #     self.low = float('-inf')
            #     self.high = 0
            elif truncation_settings == 'in-vivo':
                self.low = 0.1
                self.high = np.log(10)
            else:
                raise ValueError('`truncation_settings` ({}) not recognized'.format(
                    truncation_settings))
        elif pl.istuple(truncation_settings):
            if len(truncation_settings) != 2:
                raise ValueError('If `truncation_settings` is a tuple, it must have a ' \
                    'length of 2 ({})'.format(len(truncation_settings)))
            l,h = truncation_settings

            if (not pl.isnumeric(l)) or (not pl.isnumeric(h)):
                raise TypeError('`low` ({}) and `high` ({}) must be numerics'.format(
                    type(l), type(h)))
            if l < 0 or h < 0:
                raise ValueError('`low` ({}) and `high` ({}) must be >= 0'.format(l,h))
            if h <= l:
                raise ValueError('`low` ({}) must be strictly less than high ({})'.format(l,h))
            self.high = h
            self.low = l
        else:
            raise TypeError('`truncation_settings` ({}) type not recognized')
        self.proposal.high = self.high
        self.proposal.low = self.low

        # Set the mean
        if not pl.isstr(mean_option):
            raise TypeError('`mean_option` ({}) must be a str'.format(type(mean_option)))
        if mean_option == 'manual':
            if not pl.isnumeric(mean):
                raise TypeError('`mean` ({}) must be a numeric'.format(type(mean)))
        elif mean_option in ['auto', 'median-linear-regression']:
            # Perform linear regression
            if self.child_name == STRNAMES.GROWTH_VALUE:
                rhs = [REPRNAMES.GROWTH_VALUE, REPRNAMES.SELF_INTERACTION_VALUE]
                lhs = []
            else:
                rhs = [REPRNAMES.SELF_INTERACTION_VALUE]
                lhs = [REPRNAMES.GROWTH_VALUE]
            X = self.G.data.construct_rhs(keys=rhs,
                kwargs_dict={REPRNAMES.GROWTH_VALUE:{'with_perturbations':False}},
                index_out_perturbations=True)
            y = self.G.data.construct_lhs(keys=lhs, 
                kwargs_dict={REPRNAMES.GROWTH_VALUE:{'with_perturbations':False}},
                index_out_perturbations=True)

            prec = X.T @ X
            cov = pinv(prec, self)
            mean = cov @ X.T @ y

            if self.child_name == STRNAMES.GROWTH_VALUE:
                mean = np.median(mean[:self.G.data.n_asvs])
            else:
                mean = np.median(mean)
        else:
            raise ValueError('`mean_option` ({}) not recognized'.format(mean_option))
        self.prior.mean.override_value(mean)

        # Set the var
        if not pl.isstr(var_option):
            raise TypeError('`var_option` ({}) must be a str'.format(type(var_option)))
        if var_option == 'manual':
            if not pl.isnumeric(var):
                raise TypeError('`var` ({}) must be a numeric'.format(type(var)))
        elif var_option in ['auto', 'diffuse-linear-regression']:
            # Perform linear regression
            if self.child_name == STRNAMES.GROWTH_VALUE:
                rhs = [REPRNAMES.GROWTH_VALUE, REPRNAMES.SELF_INTERACTION_VALUE]
                lhs = []
            else:
                rhs = [REPRNAMES.SELF_INTERACTION_VALUE]
                lhs = [REPRNAMES.GROWTH_VALUE]
            X = self.G.data.construct_rhs(keys=rhs,
                kwargs_dict={REPRNAMES.GROWTH_VALUE:{'with_perturbations':False}},
                index_out_perturbations=True)
            y = self.G.data.construct_lhs(keys=lhs, 
                kwargs_dict={REPRNAMES.GROWTH_VALUE:{'with_perturbations':False}},
                index_out_perturbations=True)

            prec = X.T @ X
            cov = pinv(prec, self)
            mean = cov @ X.T @ y

            if self.child_name == STRNAMES.GROWTH_VALUE:
                mean = np.median(mean[:self.G.data.n_asvs])
            else:
                mean = np.median(mean)
            var = 1e4 * (mean**2)
        else:
            raise ValueError('`var_option` ({}) not recognized'.format(var_option))
        self.prior.var.override_value(var)

        # Set the value
        if not pl.isstr(value_option):
            raise TypeError('`value_option` ({}) must be a str'.format(type(value_option)))
        if value_option == 'manual':
            if not pl.isnumeric(value):
                raise TypeError('`value` ({}) must be a numeric'.format(type(value)))
        elif value_option in ['linear-regression']:
            # Perform linear regression
            if self.child_name == STRNAMES.GROWTH_VALUE:
                rhs = [REPRNAMES.GROWTH_VALUE, REPRNAMES.SELF_INTERACTION_VALUE]
                lhs = []
            else:
                rhs = [REPRNAMES.SELF_INTERACTION_VALUE]
                lhs = [REPRNAMES.GROWTH_VALUE]
            X = self.G.data.construct_rhs(keys=rhs,
                kwargs_dict={REPRNAMES.GROWTH_VALUE:{'with_perturbations':False}},
                index_out_perturbations=True)
            y = self.G.data.construct_lhs(keys=lhs, 
                kwargs_dict={REPRNAMES.GROWTH_VALUE:{'with_perturbations':False}},
                index_out_perturbations=True)

            prec = X.T @ X
            cov = pinv(prec, self)
            mean = cov @ X.T @ y

            if self.child_name == STRNAMES.GROWTH_VALUE:
                value = mean[:self.G.data.n_asvs]
            else:
                value = mean
        elif value_option in ['auto', 'prior-mean']:
            value = self.prior.mean.value
        else:
            raise ValueError('`value_option` ({}) not recognized'.format(value_option))
        self.value = value

        # Set the proposal variance
        if not pl.isstr(proposal_option):
            raise TypeError('`proposal_option` ({}) must be a str'.format(
                type(proposal_option)))
        elif proposal_option == 'manual':
            if not pl.isnumeric(proposal_var):
                raise TypeError('`proposal_var` ({}) must be a numeric'.format(
                    type(proposal_var)))
            if proposal_var <= 0:
                raise ValueError('`proposal_var` ({}) not proper'.format(proposal_var))
        elif proposal_option in ['auto']:
            proposal_var = (self.value ** 2)/10
        else:
            raise ValueError('`proposal_option` ({}) not recognized'.format(
                proposal_option))
        self.proposal.var.value = proposal_var

    def update_var(self):
        '''Update the `var` parameter so that we adjust the acceptance
        rate to `target_acceptance_rate`
        '''
        if self.sample_iter == 0:
            self.temp_acceptances = 0
            self.acceptances = np.zeros(self.G.inference.n_samples, dtype=bool)
        
        elif self.sample_iter > self.end_tune:
            # Don't do any more updates
            return
        
        elif self.sample_iter % self.tune == 0:
            # Update var
            acceptance_rate = self.temp_acceptances / self.tune
            if acceptance_rate > self.target_acceptance_rate:
                self.proposal.var.value *= 1.5
            else:
                self.proposal.var.value /= 1.5
            self.temp_acceptances = 0

    def update(self):
        '''First we check if we need to tune the var, which we do during
        the first half of burnin. We calculate the likelihoods in logspace
        '''
        if self.sample_iter < self.delay:
            return
        self.update_var()
        proposal_std = np.sqrt(self.proposal.var.value)

        # Get necessary data of the respective parameter
        variable = self.G[self.child_name]
        x = variable.value.ravel()
        std = np.sqrt(variable.prior.var.value)

        low = variable.low
        high = variable.high

        # propose a new value for the mean
        prev_mean = self.value
        self.proposal.mean.value = self.value
        new_mean = self.proposal.sample() # Sample a new value

        # Calculate the target distribution ll
        prev_target_ll = pl.random.truncnormal.logpdf( 
            value=prev_mean, mean=self.prior.mean.value, 
            std=np.sqrt(self.prior.var.value), low=self.low,
            high=self.high)
        for i in range(len(x)):
            prev_target_ll += pl.random.truncnormal.logpdf(
                value=x[i], mean=prev_mean, std=std,
                low=low, high=high)
        new_target_ll = pl.random.truncnormal.logpdf( 
            value=new_mean, mean=self.prior.mean.value, 
            std=np.sqrt(self.prior.var.value), low=self.low,
            high=self.high)
        for i in range(len(x)):
            new_target_ll += pl.random.truncnormal.logpdf(
                value=x[i], mean=new_mean, std=std,
                low=low, high=high)

        # Normalize by the ll of the proposal
        prev_prop_ll = pl.random.truncnormal.logpdf(
            value=prev_mean, mean=new_mean, std=proposal_std,
            low=low, high=high)
        
        new_prop_ll = pl.random.truncnormal.logpdf(
            value=new_mean, mean=prev_mean, std=proposal_std,
            low=low, high=high)

        # Accept or reject
        r = (new_target_ll - prev_prop_ll) - \
            (prev_target_ll - new_prop_ll)
        u = np.log(pl.random.misc.fast_sample_standard_uniform())

        # print('\n\n\n{} prior_mean\n----------'.format(self.child_name))
        # print('x', x)
        # print('prev_mean', prev_mean)
        # print('prev_target_ll', prev_target_ll)
        # print('prev_prop_ll', prev_prop_ll)
        # print('new mean', new_mean)
        # print('new_target_ll', new_target_ll)
        # print('new_prop_ll', new_prop_ll)
        # print('\nr', r, u)

        if r >= u:
            self.acceptances[self.sample_iter] = True
            self.value = new_mean
            self.temp_acceptances += 1
        else:
            self.value = prev_mean

    def visualize(self, path, f, section='posterior'):
        '''Render the traces in the folder `basepath` and write the 
        learned values to the file `f`.

        Parameters
        ----------
        path : str
            This is the path to write the files to
        f : _io.TextIOWrapper
            File that we are writing the values to
        section : str
            Section of the trace to compute on. Options:
                'posterior' : posterior samples
                'burnin' : burn-in samples
                'entire' : both burn-in and posterior samples

        Returns
        -------
        _io.TextIOWrapper
        '''
        f.write('\n\n###################################\n')
        f.write(self.name)
        f.write('\n###################################\n')
        if not self.G.inference.is_being_traced(self):
            f.write('`{}` not learned\n\tValue: {}\n'.format(self.name, self.value))
            return f
        summ = pl.summary(self, section=section)
        for k,v in summ.items():
            f.write('\t{}: {}\n'.format(k,v))

        # Plot the traces
        ax1, ax2 = visualization.render_trace(var=self, plt_type='both', section=section,
            include_burnin=True, log_scale=True, rasterized=True)

        # Plot the prior over the posterior
        l,h = ax1.get_xlim()
        xs = np.arange(l,h,step=(h-l)/100) 
        ys = []
        for x in xs:
            ys.append(pl.random.sics.pdf(value=x, 
                dof=self.prior.dof.value,
                scale=self.prior.scale.value))
        ax1.plot(xs, ys, label='prior', alpha=0.5, color='red', rasterized=True)
        ax1.legend()

        # Plot the acceptance rate over the trace
        ax3 = ax2.twinx()
        ax3 = visualization.render_acceptance_rate_trace(var=self, ax=ax3, 
            label='Acceptance Rate', color='red', scatter=False, rasterized=True)
        ax3.legend()
        fig = plt.gcf()
        fig.tight_layout()
        fig.suptitle(self.name)
        plt.savefig(path)
        plt.close()

        return f


class Growth(pl.variables.TruncatedNormal):
    '''Growth values of Lotka-Voltera
    '''
    def __init__(self, prior, **kwargs):
        kwargs['name'] = STRNAMES.GROWTH_VALUE
        pl.variables.TruncatedNormal.__init__(self, mean=None, var=None, low=0.,
            high=float('inf'), dtype=float, **kwargs)
        self.set_value_shape(shape=(len(self.G.data.asvs),))
        self.add_prior(prior)
        self.delay = 0
        self._initialized = False

    def __str__(self):
        return str(self.value)

    def update_str(self):
        return

    def initialize(self, value_option, truncation_settings,
        value=None, delay=0, mean=None):
        '''Initialize the growth values and hyperparamters

        Parameters
        ----------
        value_option : str
            How to initialize the values.
            Options:
                'manual'
                    Set the values manually. `value` must also be specified.
                'linear regression'
                    Set the values of the growth using linear regression
                'ones'
                    Set all of the values to 1.
                'auto'
                    Alias for 'ones'
                'prior-mean'
                    Set to the mean of the prior
        value : array
            Only necessary if `value_option` is 'manual'
        delay : int
            How many MCMC iterations to delay starting to update
        truncation_settings : str, tuple, None
            These are the settings of how you set the upper and lower limit of the
            truncated distribution. If it is None, it will default to 'standard'.
            Options
                'positive', None
                    Only constrains the values to being positive
                    low=0., high=float('inf')
                'in-vivo', 'auto'
                    Tighter constraint on the growth values.
                    low=0.1, high=ln(10)
                    These values have the following meaning
                        The slowest growing microbe will grow an order of magnitude in ~10 days
                        The fastest growing microbe will grow an order of magnitude in 1 day
                tuple(low, high)
                    These are manually specified values for the low and high
        '''
        self._initialized = True
        self._there_are_perturbations = self.G.perturbations is not None
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        self.delay = delay

        # Truncation settings
        if truncation_settings is None:
            truncation_settings = 'positive'
        if pl.isstr(truncation_settings):
            if truncation_settings == 'positive':
                self.low = 0.
                self.high = float('inf')
            elif truncation_settings in ['in-vivo', 'auto']:
                self.low = 0.1
                self.high = math.log(10)
            else:
                raise ValueError('`truncation_settings` ({}) not recognized'.format(
                    truncation_settings))
        elif pl.istuple(truncation_settings):
            if len(truncation_settings) != 2:
                raise ValueError('If `truncation_settings` is a tuple, it must have a ' \
                    'length of 2 ({})'.format(len(truncation_settings)))
            l,h = truncation_settings

            if (not pl.isnumeric(l)) or (not pl.isnumeric(h)):
                raise TypeError('`low` ({}) and `high` ({}) must be numerics'.format(
                    type(l), type(h)))
            if l < 0 or h < 0:
                raise ValueError('`low` ({}) and `high` ({}) must be >= 0'.format(l,h))
            if h <= l:
                raise ValueError('`low` ({}) must be strictly less than high ({})'.format(l,h))
            self.high = h
            self.low = l
        else:
            raise TypeError('`truncation_settings` ({}) type not recognized')

        # Setting the value
        if not pl.isstr(value_option):
            raise TypeError('`value_option` ({}) must be a str'.format(type(value_option)))
        if value_option == 'manual':
            if not pl.isarray(value):
                value = np.ones(len(self.G.data.asvs))*value
            if len(value) != self.G.data.n_asvs:
                raise ValueError('`value` ({}) must be ({}) long'.format(
                    len(value), len(self.G.data.asvs)))
            self.value = value
        elif value_option == 'linear-regression':
            rhs = [
                REPRNAMES.GROWTH_VALUE,
                REPRNAMES.SELF_INTERACTION_VALUE
            ]
            lhs = []
            X = self.G.data.construct_rhs(
                keys=rhs, kwargs_dict={REPRNAMES.GROWTH_VALUE:{'with_perturbations':False}},
                index_out_perturbations=True)
            y = self.G.data.construct_lhs(keys=lhs, index_out_perturbations=True)

            prec = X.T @ X
            cov = pinv(prec, self)
            mean = (cov @ X.transpose().dot(y)).ravel()
            self.value = np.absolute(mean[:len(self.G.data.asvs)])
        elif value_option in ['auto', 'ones']:
            self.value = np.ones(len(self.G.data.asvs), dtype=float)
        elif value_option == 'prior-mean':
            self.value = self.prior.mean.value * np.ones(self.G.data.n_asvs)
        else:
            raise ValueError('`value_option` ({}) not recognized'.format(value_option))

        logging.info('Growth value initialization: {}'.format(self.value))
        logging.info('Growth prior mean: {}'.format(self.prior.mean.value))
        logging.info('Growth truncation settings: {}'.format((self.low, self.high)))

    def update(self):
        '''Update the values using a truncated normal
        '''
        if self.sample_iter < self.delay:
            return

        self.calculate_posterior()
        self.sample()

        if not pl.isarray(self.value):
            # This will happen if there is 1 ASV
            self.value = np.array([self.value])

        if np.any(np.isnan(self.value)):
            logging.critical('mean: {}'.format(self.mean.value))
            logging.critical('var: {}'.format(self.var.value))
            logging.critical('value: {}'.format(self.value))
            raise ValueError('`Values in {} are nan: {}'.format(self.name, self.value))

        if self._there_are_perturbations:
            # If there are perturbations then we need to update their
            # matrix because the growths changed
            self.G.data.design_matrices[REPRNAMES.PERT_VALUE].update_values()

    def calculate_posterior(self):
        rhs = [REPRNAMES.GROWTH_VALUE]
        if self._there_are_perturbations:
            lhs = [
                REPRNAMES.SELF_INTERACTION_VALUE,
                REPRNAMES.CLUSTER_INTERACTION_VALUE]
        else:
            lhs = [
                REPRNAMES.SELF_INTERACTION_VALUE,
                REPRNAMES.CLUSTER_INTERACTION_VALUE]
        X = self.G.data.construct_rhs(keys=rhs,
            kwargs_dict={REPRNAMES.GROWTH_VALUE:{
                'with_perturbations':self._there_are_perturbations}})
        y = self.G.data.construct_lhs(keys=lhs)
        # X = X.toarray()

        process_prec = self.G[REPRNAMES.PROCESSVAR].build_matrix(
            cov=False, sparse=True)

        prior_prec = build_prior_covariance(G=self.G, cov=False,
            order=rhs, sparse=True)
        prior_mean = build_prior_mean(G=self.G, order=rhs).reshape(-1,1)
        pm = prior_prec @ prior_mean

        prec = X.T @ process_prec @ X + prior_prec
        cov = pinv(prec, self)

        self.mean.value = np.asarray(cov @ (X.T @ process_prec.dot(y) + pm)).ravel()
        self.var.value = np.diag(cov)

    def visualize(self, basepath, f, section='posterior', asv_formatter='%(name)s', 
        true_value=None):
        '''Render the traces in the folder `basepath` and write the 
        learned values to the file `f`.

        Parameters
        ----------
        basepath : str
            This is the loction to write the files to
        f : _io.TextIOWrapper
            File that we are writing the values to
        section : str
            Section of the trace to compute on. Options:
                'posterior' : posterior samples
                'burnin' : burn-in samples
                'entire' : both burn-in and posterior samples
        true_value : np.ndarray
            Ground truth values of the variable

        Returns
        -------
        _io.TextIOWrapper
        '''
        f.write('\n\n###################################\n')
        f.write(self.name)
        f.write('\n###################################\n')
        if not self.G.inference.is_being_traced(self):
            f.write('`{}` not learned\n\tValue: {}\n'.format(self.name, self.value))
            return f

        asvs = self.G.data.subjects.asvs
        summ = pl.summary(self, section=section)
        for key,arr in summ.items():
            f.write('{}\n'.format(key))
            for idx,ele in enumerate(arr):
                prefix = ''
                if asv_formatter is not None:
                    prefix = pl.asvname_formatter(format=asv_formatter, asv=asvs[idx], asvs=asvs)
                f.write('\t' + prefix + '{}\n'.format(ele)) 

        if section == 'posterior':
            len_posterior = self.G.inference.sample_iter + 1 - self.G.inference.burnin
        elif section == 'burnin':
            len_posterior = self.G.inference.burnin
        else:
            len_posterior = self.G.inference.sample_iter + 1

        # Plot the prior on top of the posterior
        if self.G.tracer.is_being_traced(STRNAMES.PRIOR_MEAN_GROWTH):
            prior_mean_trace = self.G[STRNAMES.PRIOR_MEAN_GROWTH].get_trace_from_disk(
                    section=section)
        else:
            prior_mean_trace = self.prior.mean.value * np.ones(len_posterior, dtype=float)
        if self.G.tracer.is_being_traced(STRNAMES.PRIOR_VAR_GROWTH):
            prior_std_trace = np.sqrt(
                self.G[STRNAMES.PRIOR_VAR_GROWTH].get_trace_from_disk(section=section))
        else:
            prior_std_trace = np.sqrt(self.prior.var.value) * np.ones(len_posterior, dtype=float)

        for idx in range(len(asvs)):
            fig = plt.figure()
            ax_posterior = fig.add_subplot(1,2,1)
            visualization.render_trace(var=self, idx=idx, plt_type='hist',
                label=section, color='blue', ax=ax_posterior, section=section,
                include_burnin=True, rasterized=True)

            # Get the limits and only look at the posterior within 20% range +- of
            # this number
            low_x, high_x = ax_posterior.get_xlim()

            arr = np.zeros(len(prior_std_trace), dtype=float)
            for i in range(len(prior_std_trace)):
                arr[i] = pl.random.truncnormal.sample(mean=prior_mean_trace[i], std=prior_std_trace[i], 
                    low=self.low, high=self.high)
            visualization.render_trace(var=arr, plt_type='hist', 
                label='prior', color='red', ax=ax_posterior, rasterized=True)

            if true_value is not None:
                ax_posterior.axvline(x=true_value[idx], color='red', alpha=0.65, 
                    label='True Value')

            ax_posterior.legend()
            ax_posterior.set_xlim(left=low_x*.8, right=high_x*1.2)

            # plot the trace
            ax_trace = fig.add_subplot(1,2,2)
            visualization.render_trace(var=self, idx=idx, plt_type='trace', 
                ax=ax_trace, section=section, include_burnin=True, rasterized=True)

            if true_value is not None:
                ax_trace.axhline(y=true_value[idx], color='red', alpha=0.65, 
                    label='True Value')
                ax_trace.legend()

            if asv_formatter is not None:
                asvname = pl.asvname_formatter(
                    format=asv_formatter,
                    asv=asvs[idx],
                    asvs=asvs)
            else:
                asvname = asvs[idx].name
            asvname = asvname.replace('/', '_').replace(' ', '_')

            fig.suptitle('Growth {}'.format(asvname))
            fig.tight_layout()
            fig.subplots_adjust(top=0.85)
            plt.savefig(os.path.join(basepath, '{}.pdf'.format(asvs[idx].name)))
            plt.close()

        return f


class SelfInteractions(pl.variables.TruncatedNormal):
    '''self-interactions of Lotka-Voltera

    Since our dynamics subtract this parameter, this parameter must be positive
    '''
    def __init__(self, prior, **kwargs):
        kwargs['name'] = STRNAMES.SELF_INTERACTION_VALUE
        pl.variables.TruncatedNormal.__init__(self, mean=None, var=None, low=0.,
            high=float('inf'), dtype=float, **kwargs)
        self.set_value_shape(shape=(len(self.G.data.asvs),))
        self.add_prior(prior)

    def __str__(self):
        return str(self.value)

    def update_str(self):
        return

    def initialize(self, value_option, truncation_settings,
        value=None, delay=0, mean=None, q=None, rescale_value=None):
        '''Initialize the self-interactions values and hyperparamters

        Parameters
        ----------
        value_option : str
            How to initialize the values.
            Options:
               'manual'
                    Set the values manually. `value` must also be specified.
                'fixed-growth'
                    Fix the growth values and then sample the self-interactions
                'strict-enforcement-partial'
                    Do an unregularized regression then take the absolute value of the numbers.
                    We assume there are no interactions and we index out the time points that have
                    perturbations in them. We assume that we do not know the growths (the growths
                    are being regressed as well).
                'strict-enforcement-full'
                    Do an unregularized regression then take the absolute value of the numbers.
                    We assume there are no interactions and we index out the time points that have
                    perturbations in them. We assume that we know the growths (the growths are on
                    the lhs)
                'steady-state', 'auto'
                    Set to the steady state values. Must also provide the quantile with the
                    parameter `q`. In here we assume that the steady state is the `q`th quantile
                    of the off perturbation data
                'prior-mean'
                    Set the value to the mean of the prior
        truncation_settings : str, 2-tuple
            How to set the truncations for the normal distribution
            (low,high)
                These are the low and high values
            'negative'
                Truncated (-inf, 0)
            'positive', 'auto'
                Truncated (0, inf)
            'human'
                This assumes that the range of the steady state abundances of the human gut
                fluctuate between 1e2 and  1e13. We requie that the growth value be initialized first.
                We set the vlaues to be (growth.high/1e14, growth.low/1e2)
            'mouse'
                This assumes that the range of the steady state abundances of the mouse gut
                fluctuate between 1e2 and  1e12. We requie that the growth value be initialized first.
                We set the vlaues to be (growth.high/1e13, growth.low/1e2)
        value : array
            Only necessary if `value_option` is 'manual'
        mean : array
            Only necessary if `mean_option` is 'manual'
        delay : int
            How many MCMC iterations to delay starting to update
        rescale_value : None, float
            This is the rescale value of the qPCR. This will rescale the truncation settings.
            This is only used for either the 'mouse' or 'human' settings
        '''
        self._there_are_perturbations = self.G.perturbations is not None
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        self.delay = delay

        # Set truncation settings
        if pl.isstr(truncation_settings):
            # if truncation_settings == 'negative':
            #     self.low = float('-inf')
            #     self.high= 0
            if truncation_settings in ['mouse', 'human']:
                growth = self.G[STRNAMES.GROWTH_VALUE]
                if not growth._initialized:
                    raise ValueError('Growth values `{}` must be initialized first'.format(
                        STRNAMES.GROWTH_VALUE))
                if truncation_settings == 'mouse':
                    high = 1e13
                else:
                    high = 1e14
                low = 1e2
                if rescale_value is not None:
                    if not pl.isnumeric(rescale_value):
                        raise TypeError('`rescale_value` ({}) must be a numeric'.format(
                            type(rescale_value)))
                    if rescale_value <= 0:
                        raise ValueError('`rescale_value` ({}) must be > 0'.format(rescale_value))
                    high *= rescale_value
                    low *= rescale_value
                self.low = growth.high/low
                self.high = growth.low/high
            elif truncation_settings in ['auto', 'positive']:
                self.low = 0
                self.high = float('inf')
            else:
                raise ValueError('`truncation_settings) ({}) not recognized'.format(
                    truncation_settings))
        elif pl.istuple(truncation_settings):
            if len(truncation_settings) != 2:
                raise ValueError('If `truncation_settings` is a tuple, it must have a ' \
                    'length of 2 ({})'.format(len(truncation_settings)))
            l,h = truncation_settings

            if (not pl.isnumeric(l)) or (not pl.isnumeric(h)):
                raise TypeError('`low` ({}) and `high` ({}) must be numerics'.format(
                    type(l), type(h)))
            if l < 0 or h < 0:
                raise ValueError('`low` ({}) and `high` ({}) must be >= 0'.format(l,h))
            if h <= l:
                raise ValueError('`low` ({}) must be strictly less than high ({})'.format(l,h))
            self.high = h
            self.low = l
        else:
            raise TypeError('`truncation_settings` ({}) must be a tuple or str'.format(
                type(truncation_settings)))

        # Set value option
        if not pl.isstr(value_option):
            raise TypeError('`value_option` ({}) must be a str'.format(type(value_option)))
        if value_option == 'manual':
            if not pl.isarray(value):
                value = np.ones(len(self.G.data.asvs))*value
            if len(value) != self.G.data.n_asvs:
                raise ValueError('`value` ({}) must be ({}) long'.format(
                    len(value), len(self.G.data.asvs)))
            self.value = value
        elif value_option == 'fixed-growth':
            X = self.G.data.construct_rhs(keys=[REPRNAMES.SELF_INTERACTION_VALUE],
                index_out_perturbations=True)
            y = self.G.data.construct_lhs(keys=[REPRNAMES.GROWTH_VALUE], kwargs_dict={
                REPRNAMES.GROWTH_VALUE:{'with_perturbations':False}},
                index_out_perturbations=True)
            prec = X.T @ X
            cov = pinv(prec, self)
            self.value = np.absolute((cov @ X.transpose().dot(y)).ravel())
        elif 'strict-enforcement' in value_option:
            if 'full' in value_option:
                rhs = [REPRNAMES.SELF_INTERACTION_VALUE]
                lhs = [REPRNAMES.GROWTH_VALUE]
            elif 'partial' in value_option:
                lhs = []
                rhs = [REPRNAMES.GROWTH_VALUE, REPRNAMES.SELF_INTERACTION_VALUE]
            else:
                raise ValueError('`value_option` ({}) not recognized'.format(value_option))
            X = self.G.data.construct_rhs(
                keys=rhs, kwargs_dict={REPRNAMES.GROWTH_VALUE:{'with_perturbations':False}},
                index_out_perturbations=True)
            y = self.G.data.construct_lhs(keys=lhs, index_out_perturbations=True)

            prec = X.T @ X
            cov = pinv(prec, self)
            mean = (cov @ X.transpose().dot(y)).ravel()
            self.value = np.absolute(mean[len(self.G.data.asvs):])
        elif value_option == 'prior-mean':
            self.value = self.prior.mean.value * np.ones(self.G.data.n_asvs)
        elif value_option in ['steady-state', 'auto']:
            # check quantile
            if not pl.isnumeric(q):
                raise TypeError('`q` ({}) must be numeric'.format(type(q)))
            if q < 0 or q > 1:
                raise ValueError('`q` ({}) must be [0,1]'.format(q))

            # Get the data off perturbation
            datas = None
            for ridx in range(self.G.data.n_replicates):
                if self._there_are_perturbations:
                    # Exclude the data thats in a perturbation
                    base_idx = 0
                    for start,end in self.G.data.tidxs_in_perturbation[ridx]:
                        if datas is None:
                            datas = self.G.data.data[ridx][:,base_idx:start]
                        else:
                            datas = np.hstack((datas, self.G.data.data[ridx][:,base_idx:start]))
                        base_idx = end
                    if end != self.G.data.data[ridx].shape[1]:
                        datas = np.hstack((datas, self.G.data.data[ridx][:,base_idx:]))
                else:
                    if datas is None:
                        datas = self.G.data.data[ridx]
                    else:
                        datas = np.hstack((datas, self.G.data.data[ridx]))

            # Set the steady-state for each ASV
            ss = np.quantile(datas, q=q, axis=1)

            # Get the self-interactions by using the values of the growth terms
            self.value = 1/ss
        elif value_option == 'linear-regression':
            
            rhs = [REPRNAMES.SELF_INTERACTION_VALUE]
            lhs = [REPRNAMES.GROWTH_VALUE]
            X = self.G.data.construct_rhs(keys=rhs,
                index_out_perturbations=True)
            y = self.G.data.construct_lhs(keys=lhs, index_out_perturbations=True,
                kwargs_dict={REPRNAMES.GROWTH_VALUE:{'with_perturbations':False}})

            prec = X.T @ X
            cov = pinv(prec, self)
            mean = cov @ X.T @ y
            self.value = np.asarray(mean).ravel()
        else:
            raise ValueError('`value_option` ({}) not recognized'.format(value_option))

        logging.info('Self-interactions value initialization: {}'.format(self.value))
        logging.info('Self-interactions truncation settings: {}'.format((self.low, self.high)))

    def update(self):
        if self.sample_iter < self.delay:
            return

        self.calculate_posterior()
        self.sample()

        if not pl.isarray(self.value):
            # This will happen if there is 1 ASV
            self.value = np.array([self.value])

        if np.any(np.isnan(self.value)):
            logging.critical('mean: {}'.format(self.mean.value))
            logging.critical('var: {}'.format(self.var.value))
            logging.critical('value: {}'.format(self.value))
            raise ValueError('`Values in {} are nan: {}'.format(self.name, self.value))

    def calculate_posterior(self):

        rhs = [REPRNAMES.SELF_INTERACTION_VALUE]
        if self._there_are_perturbations:
            lhs = [
                REPRNAMES.GROWTH_VALUE,
                REPRNAMES.CLUSTER_INTERACTION_VALUE]
        else:
            lhs = [
                REPRNAMES.GROWTH_VALUE,
                REPRNAMES.CLUSTER_INTERACTION_VALUE]
        X = self.G.data.construct_rhs(keys=rhs)
        y = self.G.data.construct_lhs(keys=lhs, kwargs_dict={REPRNAMES.GROWTH_VALUE:{
                'with_perturbations':self._there_are_perturbations}})
        process_prec = self.G[REPRNAMES.PROCESSVAR].build_matrix(
            cov=False, sparse=True)
        prior_prec = build_prior_covariance(G=self.G, cov=False,
            order=rhs, sparse=True)

        pm = prior_prec @ (self.prior.mean.value * np.ones(self.G.data.n_asvs).reshape(-1,1))

        prec = X.T @ process_prec @ X + prior_prec
        cov = pinv(prec, self)
        self.mean.value = np.asarray(cov @ (X.T @ process_prec.dot(y) + pm)).ravel()
        self.var.value = np.diag(cov)

    def visualize(self, basepath, f, section='posterior', 
        asv_formatter=pl.ASVNAME_PAPER_FORMAT, true_value=None):
        '''Render the traces in the folder `basepath` and write the 
        learned values to the file `f`.

        Parameters
        ----------
        basepath : str
            This is the loction to write the files to
        f : _io.TextIOWrapper
            File that we are writing the values to
        section : str
            Section of the trace to compute on. Options:
                'posterior' : posterior samples
                'burnin' : burn-in samples
                'entire' : both burn-in and posterior samples
        true_value : np.ndarray
            Ground truth values of the variable

        Returns
        -------
        _io.TextIOWrapper
        '''
        f.write('\n\n###################################\n')
        f.write(self.name)
        f.write('\n###################################\n')
        if not self.G.inference.is_being_traced(self):
            f.write('`{}` not learned\n\tValue: {}\n'.format(self.name, self.value))
            return f

        asvs = self.G.data.subjects.asvs
        summ = pl.summary(self, section=section)
        for key,arr in summ.items():
            f.write('{}\n'.format(key))
            for idx,ele in enumerate(arr):
                prefix = ''
                if asv_formatter is not None:
                    prefix = pl.asvname_formatter(format=asv_formatter, asv=asvs[idx], asvs=asvs)
                f.write('\t' + prefix + '{}\n'.format(ele)) 

        if section == 'posterior':
            len_posterior = self.G.inference.sample_iter + 1 - self.G.inference.burnin
        elif section == 'burnin':
            len_posterior = self.G.inference.burnin
        else:
            len_posterior = self.G.inference.sample_iter + 1

        # Plot the prior on top of the posterior
        if self.G.tracer.is_being_traced(STRNAMES.PRIOR_MEAN_SELF_INTERACTIONS):
            prior_mean_trace = self.G[STRNAMES.PRIOR_MEAN_SELF_INTERACTIONS].get_trace_from_disk(
                    section=section)
        else:
            prior_mean_trace = self.prior.mean.value * np.ones(len_posterior, dtype=float)
        if self.G.tracer.is_being_traced(STRNAMES.PRIOR_VAR_SELF_INTERACTIONS):
            prior_std_trace = np.sqrt(
                self.G[STRNAMES.PRIOR_VAR_SELF_INTERACTIONS].get_trace_from_disk(section=section))
        else:
            prior_std_trace = np.sqrt(self.prior.var.value) * np.ones(len_posterior, dtype=float)

        for idx in range(len(asvs)):
            fig = plt.figure()
            ax_posterior = fig.add_subplot(1,2,1)
            visualization.render_trace(var=self, idx=idx, plt_type='hist',
                label=section, color='blue', ax=ax_posterior, section=section,
                include_burnin=True, rasterized=True, log_scale=True)

            # Get the limits and only look at the posterior within 20% range +- of
            # this number
            low_x, high_x = ax_posterior.get_xlim()

            arr = np.zeros(len(prior_std_trace), dtype=float)
            for i in range(len(prior_std_trace)):
                arr[i] = pl.random.truncnormal.sample(mean=prior_mean_trace[i], std=prior_std_trace[i], 
                    low=self.low, high=self.high)
            visualization.render_trace(var=arr, plt_type='hist', log_scale=True,
                label='prior', color='red', ax=ax_posterior, rasterized=True)

            if true_value is not None:
                ax_posterior.axvline(x=true_value[idx], color='red', alpha=0.65, 
                    label='True Value')

            ax_posterior.legend()
            ax_posterior.set_xlim(left=low_x*.8, right=high_x*1.2)

            # plot the trace
            ax_trace = fig.add_subplot(1,2,2)
            visualization.render_trace(var=self, idx=idx, plt_type='trace', 
                ax=ax_trace, section=section, include_burnin=True, rasterized=True,
                log_scale=True)

            if true_value is not None:
                ax_trace.axhline(y=true_value[idx], color='red', alpha=0.65, 
                    label='True Value')
                ax_trace.legend()

            if asv_formatter is not None:
                asvname = pl.asvname_formatter(
                    format=asv_formatter,
                    asv=asvs[idx],
                    asvs=asvs)
            else:
                asvname = asvs[idx].name
            asvname = asvname.replace('/', '_').replace(' ', '_')

            fig.suptitle('Self-Interactions {}'.format(asvname))
            fig.tight_layout()
            fig.subplots_adjust(top=0.85)
            plt.savefig(os.path.join(basepath, '{}.pdf'.format(asvs[idx].name)))
            plt.close()

        return f
    

class RegressCoeff(pl.variables.MVN):
    '''This is the posterior of the regression coefficients.
    The current posterior assumes a prior mean of 0.

    This class samples the growth, self-interactions, and cluster
    interactions jointly.

    Parameters
    ----------
    growth : posterior.Growth
        This is the class that has the growth variables
    self_interactions : posterior.SelfInteractions
        The self interaction terms for the ASVs
    interactions : ClusterInteractionValue
        These are the cluster interaction values
    pert_mag : PerturbationMagnitudes, None
        These are the magnitudes of the perturbation parameters (per clsuter)
        Set to None if there are no perturbations
    '''
    def __init__(self, growth, self_interactions, interactions,
        pert_mag, **kwargs):

        if not issubclass(growth.__class__, Growth):
            raise ValueError('`growth` ({}) must be a subclass of the Growth ' \
                'class'.format(type(growth)))
        if not issubclass(self_interactions.__class__, SelfInteractions):
            raise ValueError('`self_interactions` ({}) must be a subclass of the SelfInteractions ' \
                'class'.format(type(self_interactions)))
        if not issubclass(interactions.__class__, ClusterInteractionValue):
            raise ValueError('`interactions` ({}) must be a subclass of the Interactions ' \
                'class'.format(type(interactions)))
        if pert_mag is not None:
            if not issubclass(pert_mag.__class__, PerturbationMagnitudes):
                raise ValueError('`pert_mag` ({}) must be a subclass of the PerturbationMagnitudes ' \
                    'class'.format(type(pert_mag)))


        kwargs['name'] = STRNAMES.REGRESSCOEFF
        pl.variables.MVN.__init__(self, mean=None, cov=None, dtype=float, **kwargs)

        self.n_asvs = self.G.data.n_asvs
        self.growth = growth
        self.self_interactions = self_interactions
        self.interactions = interactions
        self.pert_mag = pert_mag
        self.clustering = interactions.clustering

        # These serve no functional purpose but we do it so that they are
        # connected in the graph structure. Each of these should have their
        # prior already initialized
        self.add_parent(self.growth)
        self.add_parent(self.self_interactions)
        self.add_parent(self.interactions)

    def __str__(self):
        '''Make it more readable
        '''
        try:
            a = 'Growth:\n{}\nSelf Interactions:\n{}\nInteractions:\n{}\nPerturbations:\n{}\n' \
                'Acceptances:\n{}'.format(
                self.growth.value, self.self_interactions.value,
                str(self.G[REPRNAMES.CLUSTER_INTERACTION_VALUE]),
                str(self.pert_mag), np.mean(
                    self.acceptances[ np.max([self.sample_iter-50, 0]):self.sample_iter], axis=0))
        except:
            a = 'Growth:\n{}\nSelf Interactions:\n{}\nInteractions:\n{}\nPerturbations:\n{}'.format(
                self.growth.value, self.self_interactions.value,
                str(self.G[REPRNAMES.CLUSTER_INTERACTION_VALUE]),
                str(self.pert_mag))
        return a

    def initialize(self, update_jointly_pert_inter, update_jointly_growth_si, 
        tune=None, end_tune=None):
        '''The interior objects are initialized by themselves. Define which variables
        get updated together.

        Note that the interactions and perturbations will always be updated before the
        growth rates and the self-interactions

        Interactions and perturbations
        ------------------------------
        These are conjugate and have a normal prior. If these are said to be updated
        jointly then we can sample directly with Gibbs sampling.

        Growths and self-interactions
        -----------------------------
        These are conjugate and have a truncated normal prior. If they are set to be 
        updated together, then we must do MH because we cannot sample from a truncated
        multivariate gaussian.

        Parameters
        ----------
        update_jointly_pert_inter : bool
            If True, update the interactions and the perturbations jointly.
            If False, update the interactions and perturbations separately - you
            randomly choose which one to update first.
        update_jointly_growth_si : bool
            If True, update the interactions and the perturbations jointly.
            If False, update the interactions and perturbations separately - you
            randomly choose which one to update first.
        '''
        self._there_are_perturbations = self.G.perturbations is not None
        if not pl.isbool(update_jointly_growth_si):
            raise TypeError('`update_jointly_growth_si` ({}) must be a bool'.format(
                type(update_jointly_growth_si)))
        if not pl.isbool(update_jointly_pert_inter):
            raise TypeError('`update_jointly_pert_inter` ({}) must be a bool'.format(
                type(update_jointly_pert_inter)))

        self.update_jointly_growth_si = update_jointly_growth_si
        self.update_jointly_pert_inter = update_jointly_pert_inter
        self.sample_iter = 0

        if self.update_jointly_growth_si:
            raise NotImplementedError('Not Implemented')
                
    # @profile
    def asarray(self):
        '''
        Builds the full regression coefficient vector. If `asv_id` and
        `cid` are None, build the entire thing. Else build it for
        the ASV or cluster specifically.

        Parameters
        ----------
        '''
        # build the entire thing
        a = np.append(self.growth.value, self.self_interactions.value)
        a = np.append(a, self.interactions.obj.get_values(use_indicators=True))
        return a

    def update(self):
        '''Either updated jointly using multivariate normal or update independently
        using truncated normal distributions for growth and self-interactions.

        Always update the one that the interactions is in first
        '''
        self._update_perts_and_inter()
        if self._there_are_perturbations:
            self.G.data.design_matrices[REPRNAMES.GROWTH_VALUE].build_with_perturbations()
        
        self._update_growth_and_self_interactions()
        self.sample_iter += 1

        if self._there_are_perturbations:
            # If there are perturbations then we need to update their
            # matrix because the growths changed
            self.G.data.design_matrices[REPRNAMES.PERT_VALUE].update_values()

    # @profile
    def _update_perts_and_inter(self):
        '''Update the with Gibbs sampling of a multivariate normal.

        Parameters
        ----------
        args : tuple
            This is a tuple of length > 1 that holds the variables on what to update
            together
        '''
        if not self.update_jointly_pert_inter:
            # Update separately
            if pl.random.misc.fast_sample_standard_uniform() < 0.5:
                self.G[REPRNAMES.CLUSTER_INTERACTION_VALUE].update()
                if self._there_are_perturbations:
                    self.G[REPRNAMES.PERT_VALUE].update()
            else:
                if self._there_are_perturbations:
                    self.G[REPRNAMES.PERT_VALUE].update()
                self.G[REPRNAMES.CLUSTER_INTERACTION_VALUE].update()
        else:
            # Update jointly
            rhs = []
            lhs = []
            if self.interactions.obj.sample_iter >= \
                self.interactions.delay:
                rhs.append(REPRNAMES.CLUSTER_INTERACTION_VALUE)
            else:
                lhs.append(REPRNAMES.CLUSTER_INTERACTION_VALUE)
            if self._there_are_perturbations:
                if self.pert_mag.sample_iter >= self.pert_mag.delay:
                    rhs.append(REPRNAMES.PERT_VALUE)
                else:
                    lhs.append(REPRNAMES.PERT_VALUE)

            if len(rhs) == 0:
                return

            lhs += [REPRNAMES.GROWTH_VALUE, REPRNAMES.SELF_INTERACTION_VALUE]
            X = self.G.data.construct_rhs(keys=rhs)
            if X.shape[1] == 0:
                logging.info('No columns, skipping')
                return
            y = self.G.data.construct_lhs(keys=lhs,
                kwargs_dict={REPRNAMES.GROWTH_VALUE:{'with_perturbations': False}})

            process_prec = self.G[REPRNAMES.PROCESSVAR].build_matrix(
                cov=False, sparse=True)
            prior_prec = build_prior_covariance(G=self.G, cov=False,
                order=rhs, sparse=True)
            prior_means = build_prior_mean(G=self.G,order=rhs).reshape(-1,1)

            # Make the prior covariance matrix and process varaince
            prec = X.T @ process_prec @ X + prior_prec
            self.cov.value = pinv(prec, self)
            self.mean.value = np.asarray(self.cov.value @ (X.T @ process_prec.dot(y) + \
                prior_prec @ prior_means)).ravel()

            # sample posterior jointly and then assign the values to each coefficient
            # type, respectfully
            try:
                value = self.sample()
            except:
                logging.critical('failed here, updating separately')
                self.pert_mag.update()
                self.interactions.update()
                return

            i = 0
            if REPRNAMES.CLUSTER_INTERACTION_VALUE in rhs:
                l = self.interactions.obj.num_pos_indicators()
                self.interactions.value = value[:l]
                self.interactions.set_values(arr=value[:l], use_indicators=True)
                self.interactions.update_str()
                i += l
            if self._there_are_perturbations:
                if REPRNAMES.PERT_VALUE in rhs:
                    self.pert_mag.value = value[i:]
                    self.pert_mag.set_values(arr=value[i:], use_indicators=True)
                    self.pert_mag.update_str()
                    self.G.data.design_matrices[REPRNAMES.GROWTH_VALUE].update_value()
                    # self.G.data.design_matrices[REPRNAMES.PERT_VALUE].build()

    def _update_acceptances(self):
        if self.growth.sample_iter == 0:
            self.temp_acceptances= np.zeros(len(self.G.data.asvs), dtype=int)
            self.acceptances = np.zeros(shape=(self.G.inference.n_samples, 
                len(self.G.data.asvs)), dtype=bool)
        elif self.growth.sample_iter > self.end_tune:
            return
        elif self.growth.sample_iter % self.tune == 0:
            self.temp_acceptances = np.zeros(len(self.G.data.asvs), dtype=int)

    def _update_growth_and_self_interactions(self):
        '''Update the growth and self-interactions
        Our proposal is the posterior distribution.
        '''
        if not self.update_jointly_growth_si:
            # Update separately
            self.growth.update()
            self.self_interactions.update()
        else:
            # Update together
            raise NotImplementedError('Not Implemented')

    def add_trace(self):
        '''Trace values for growth, self-interactions, and cluster interaction values
        '''
        self.growth.add_trace()
        self.self_interactions.add_trace()
        self.interactions.add_trace()
        if self._there_are_perturbations:
            self.pert_mag.add_trace()

    def set_trace(self):
        self.growth.set_trace()
        self.self_interactions.set_trace()
        self.interactions.set_trace()
        if self._there_are_perturbations:
            self.pert_mag.set_trace()


# Perturbations
# -------------
class PerturbationMagnitudes(pl.variables.Normal):
    '''These update the perturbation values jointly.
    '''
    def __init__(self, **kwargs):
        '''Parameters
        '''

        kwargs['name'] = STRNAMES.PERT_VALUE
        pl.variables.Normal.__init__(self, mean=None, var=None, dtype=float, **kwargs)
        self.perturbations = self.G.perturbations

    def __str__(self):
        s = 'Perturbation Magnitudes (multiplicative)'
        for perturbation in self.perturbations:
            s += '\n\t perturbation {}: {}'.format(
                perturbation.name, perturbation.cluster_array(only_pos_ind=True))
        return s

    def __len__(self):
        '''Return the number of on indicators
        '''
        n = 0
        for perturbation in self.perturbations:
            n += perturbation.indicator.num_on_clusters()
        return n

    def set_values(self, arr, use_indicators=True):
        '''Set the values of the perturbation of them stacked one on top of each other

        Parameters
        ----------
        arr : np.ndarray
            Values for all of the perturbations in order
        use_indicators : bool
            If True, the values only refer to the on indicators
        '''
        i = 0
        for perturbation in self.perturbations:
            l = perturbation.indicator.num_on_clusters()
            perturbation.set_values_from_array(values=arr[i:i+l],
                use_indicators=use_indicators)
            i += l

    def update_str(self):
        return

    @property
    def sample_iter(self):
        return self.perturbations[0].sample_iter

    def initialize(self, value_option, value=None, mean=None, var=None, delay=0):
        '''Initialize the prior and the value of the perturbation. We assume that
        each perturbation has the same hyperparameters for the prior

        Parameters
        ----------
        value_option : str
            How to initialize the values. Options:
                'manual'
                    Set the value manually, `value` must also be specified
                'zero'
                    Set all the values to zero.
                'auto', 'prior-mean'
                    Initialize to the same value as the prior mean
        delay : int, None
            How many MCMC iterations to delay the update of the values.
        mean, var, value : int, float, array
            - Only necessary if any of the options are 'manual'
        '''
        if delay is None:
            delay = 0
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        self.delay = delay
        for perturbation in self.perturbations:
            perturbation.magnitude.set_signal_when_clusters_change(True)

        # Set the value of the perturbations
        if not pl.isstr(value_option):
            raise TypeError('`value_option` ({}) must be a str'.format(type(value_option)))
        if value_option == 'manual':
            for pidx, perturbation in enumerate(self.perturbations):
                v = value[pidx]
                for cidx, val in enumerate(v):
                    cid = perturbation.clustering.order[cidx]
                    perturbation.indicator.value[cid] = not np.isnan(val)
                    perturbation.magnitude.value[cid] = val if not np.isnan(val) else 0

        elif value_option == 'zero':
            for perturbation in self.perturbations:
                for cid in perturbation.clustering.order:
                    perturbation.magnitude.value[cid] = 0
        elif value_option in ['auto', 'prior-mean']:
            for perturbation in self.perturbations:
                mean = perturbation.magnitude.prior.mean.value
                for cid in perturbation.clustering.order:
                    perturbation.magnitude.value[cid] = mean
        else:
            raise ValueError('`value_option` ({}) not recognized'.format(value_option))


        s = 'Perturbation magnitude initialization results:\n'
        for perturbation in self.perturbations:
            if perturbation.name is not None:
                a = perturbation.name
            s += '\tPerturbation {}:\n' \
                '\t\tvalue: {}\n'.format(a, perturbation.magnitude.cluster_array())
        logging.info(s)

    def update(self):
        '''Update with a gibbs step jointly
        '''
        if self.sample_iter < self.delay:
            return

        n_on = [perturbation.indicator.num_on_clusters() for perturbation in \
            self.perturbations]
        
        if n_on == 0:
            return

        rhs = [REPRNAMES.PERT_VALUE]
        lhs = [
            REPRNAMES.GROWTH_VALUE,
            REPRNAMES.SELF_INTERACTION_VALUE,
            REPRNAMES.CLUSTER_INTERACTION_VALUE]
        X = self.G.data.construct_rhs(keys=rhs, toarray=True)
        y = self.G.data.construct_lhs(keys=lhs,
            kwargs_dict={REPRNAMES.GROWTH_VALUE:{
                'with_perturbations':False}})

        process_prec = self.G[REPRNAMES.PROCESSVAR].prec
        prior_prec = build_prior_covariance(G=self.G, cov=False, order=rhs, sparse=False)

        prior_mean = build_prior_mean(G=self.G, order=rhs).reshape(-1,1)

        a = X.T * process_prec
        prec = a @ X + prior_prec
        cov = pinv(prec, self)
        mean = np.asarray(cov @ (a @ y + prior_prec @ prior_mean)).ravel()

        # print('\n\ny\n',np.hstack((y, self.G.data.lhs.vector.reshape(-1,1))))
        # print(self.G[REPRNAMES.CLUSTER_INTERACTION_VALUE].value)
        # print(self.G[REPRNAMES.GROWTH_VALUE].value)
        # print(self.G[REPRNAMES.SELF_INTERACTION_VALUE].value)

        self.mean.value = mean
        self.var.value = np.diag(cov)
        value = self.sample()

        if np.any(np.isnan(value)):
            logging.critical('mean: {}'.format(self.mean.value))
            logging.critical('var: {}'.format(self.var.value))
            logging.critical('value: {}'.format(self.value))
            logging.critical('prior mean: {}'.format(prior_mean.ravel()))
            raise ValueError('`Values in {} are nan: {}'.format(self.name, self.value))

        i = 0
        for pidx, perturbation in enumerate(self.perturbations):
            perturbation.set_values_from_array(value[i:i+n_on[pidx]], use_indicators=True)
            i += n_on[pidx]

        # Rebuild the design matrix
        self.G.data.design_matrices[REPRNAMES.GROWTH_VALUE].build_with_perturbations()

    def set_trace(self, *args, **kwargs):
        for perturbation in self.perturbations:
            perturbation.set_trace(*args, **kwargs)

    def add_trace(self, *args, **kwargs):
        for perturbation in self.perturbations:
            perturbation.add_trace(*args, **kwargs)

    def asarray(self):
        '''Get an array of the perturbation magnitudes
        '''
        a = []
        for perturbation in self.perturbations:
            a.append(perturbation.cluster_array(only_pos_ind=True))
        return np.asarray(list(itertools.chain.from_iterable(a)))

    def toarray(self):
        return self.asarray()


class PerturbationProbabilities(pl.Node):
    '''This is the probability for a positive interaction for a perturbation
    '''
    def __init__(self, **kwargs):
        '''Parameters

        prior (pl.variables.Beta)
            - prior probability
        pert_n (int)
            - This is the perturbation number that it corresponds to
        **kwargs
            - Other options like graph, value
        '''
        kwargs['name'] = STRNAMES.PERT_INDICATOR_PROB
        pl.Node.__init__(self, **kwargs)
        self.perturbations = self.G.perturbations

    def __str__(self):
        s = 'Perturbation Indicator probabilities'
        for perturbation in self.perturbations:
            s += '\n\tperturbation {}: {}'.format(
                perturbation.name,
                perturbation.probability.value)
        return s

    @property
    def sample_iter(self):
        return self.perturbations[0].probability.sample_iter

    def initialize(self, value_option, hyperparam_option, a=None, b=None, value=None,
        N='auto', delay=0):
        '''Initialize the hyperparameters of the prior and the value. Each
        perturbation has the same prior.

        Parameters
        ----------
        value_option : str
            How to initialize the values. Options:
                'manual'
                    Set the value manually, `value` must also be specified
                'auto', 'prior-mean'
                    Initialize the value as the prior mean
        hyperparam_option : str
            How to initialize `a` and `b`. Options:
                'manual'
                    Set the value manually. `a` and `b` must also be specified
                'weak-agnostic' or 'auto'
                    a=b=0.5
                'strong-dense'
                    a = N, N are the expected number of clusters
                    b = 0.5
                'strong-sparse'
                    a = 0.5
                    b = N, N are the expected number of clusters
                'very-strong-sparse'
                    a = 0.5
                    b = N, N are the expected number of ASVs
        N : str, int
            This is the number of clusters to set the hyperparam options to 
            (if they are dependent on the number of cluster). If 'auto', set to the expected number
            of clusters from a dirichlet process. Else use this number (must be an int).
        delay : int
            How many MCMC iterations to delay starting the update of the variable
        value, a, b : int, float
            User specified values
            Only necessary if `hyperparam_option` == 'manual'
        '''
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        self.delay = delay

        # Set the hyper-parameters
        if not pl.isstr(hyperparam_option):
            raise ValueError('`hyperparam_option` ({}) must be a str'.format(type(hyperparam_option)))
        if hyperparam_option == 'manual':
            if (not pl.isnumeric(a)) or (not pl.isnumeric(b)):
                raise TypeError('If `hyperparam_option` is "manual" then `a` ({})' \
                    ' and `b` ({}) must be numerics'.format(type(a), type(b)))
        elif hyperparam_option in ['auto', 'weak-agnostic']:
            a = 0.5
            b = 0.5
        elif hyperparam_option == 'strong-dense':
            if pl.isstr(N):
                if N == 'auto':
                    N = expected_n_clusters(G=self.G)
                else:
                    raise ValueError('`N` ({}) nto recognized'.format(N))
            elif pl.isint(N):
                if N < 0:
                    raise ValueError('`N` ({}) must be positive'.format(N))
            else:
                raise TypeError('`N` ({}) type not recognized'.format(type(N)))
            a = N
            b = 0.5
        elif hyperparam_option == 'strong-sparse':
            if pl.isstr(N):
                if N == 'auto':
                    N = expected_n_clusters(G=self.G)
                else:
                    raise ValueError('`N` ({}) nto recognized'.format(N))
            elif pl.isint(N):
                if N < 0:
                    raise ValueError('`N` ({}) must be positive'.format(N))
            else:
                raise TypeError('`N` ({}) type not recognized'.format(type(N)))
            a = 0.5
            b = N
        elif hyperparam_option == 'very-strong-sparse':
            N = self.G.data.n_asvs
            a = 0.5
            b = N
        else:
            raise ValueError('`hyperparam_option` ({}) not recognized'.format(hyperparam_option))
        for perturbation in self.perturbations:
            perturbation.probability.prior.a.override_value(a)
            perturbation.probability.prior.b.override_value(b)

        # Set the value
        if not pl.isstr(value_option):
            raise TypeError('`value_option` ({}) must be a str'.format(type(value_option)))
        if value_option == 'manual':
            if not pl.isnumeric(value):
                raise TypeError('If `value_option` is "manual" then `value` ({})' \
                    ' must be a numeric'.format(type(value)))
            for perturbation in self.perturbations:
                perturbation.probability.value = value
        elif value_option in ['auto', 'prior-mean']:
            for perturbation in self.perturbations:
                perturbation.probability.value = perturbation.probability.prior.mean()
        else:
            raise ValueError('`value_option` ({}) not recognized'.format(value_option))

        s = 'Perturbation indicator probability initialization results:\n'
        for i, perturbation in enumerate(self.perturbations):
            s += '\tPerturbation {}:\n' \
                '\t\tprior a: {}\n' \
                '\t\tprior b: {}\n' \
                '\t\tvalue: {}\n'.format(i,
                    perturbation.probability.prior.a.value,
                    perturbation.probability.prior.b.value,
                    perturbation.probability.value)
        logging.info(s)

    def update(self):
        '''Update according to how many positive and negative indicators there
        are
        '''
        if self.sample_iter < self.delay:
            return
        for perturbation in self.perturbations:
            num_pos = perturbation.indicator.num_on_clusters()
            num_neg = len(perturbation.clustering.clusters) - num_pos
            perturbation.probability.a.value = perturbation.probability.prior.a.value + num_pos
            perturbation.probability.b.value = perturbation.probability.prior.b.value + num_neg
            perturbation.probability.sample()

    def set_trace(self, *args, **kwargs):
        for perturbation in self.perturbations:
            perturbation.probability.set_trace(*args, **kwargs)

    def add_trace(self, *args, **kwargs):
        for perturbation in self.perturbations:
            perturbation.probability.add_trace(*args, **kwargs)


class PerturbationIndicators(pl.Node):
    '''This is the indicator for a perturbation

    We only need to trace once for the perturbations. Our default is to only
    trace from the magnitudes. Thus, we only trace the indicators (here) if
    we are learning here and not learning the magnitudes.
    '''
    def __init__(self, need_to_trace, relative, **kwargs):
        '''Parameters
        '''
        kwargs['name'] = STRNAMES.PERT_INDICATOR
        pl.Node.__init__(self, **kwargs)
        self.need_to_trace = need_to_trace
        self.perturbations = self.G.perturbations
        self.clustering = None
        self._time_taken = None
        if relative:
            self.update = self.update_relative
        else:
            self.update = self.update_slow

    def __str__(self):
        s = 'Perturbation Indicators - time: {}s'.format(self._time_taken)
        for perturbation in self.perturbations:
            arr = perturbation.indicator.cluster_bool_array()
            s += '\nperturbation {} ({}/{}): {}'.format(perturbation.name,
                np.sum(arr), len(arr), arr)
        return s

    @property
    def sample_iter(self):
        return self.perturbations[0].sample_iter

    def add_trace(self):
        '''Only trace if perturbation indicators are being learned and the
        perturbation value is not being learned
        '''
        if self.need_to_trace:
            for perturbation in self.perturbations:
                perturbation.add_trace()

    def set_trace(self, *args, **kwargs):
        '''Only trace if perturbation indicators are being learned and the
        perturbation value is not being learned
        '''
        if self.need_to_trace:
            for perturbation in self.perturbations:
                perturbation.set_trace(*args, **kwargs)

    def initialize(self, value_option, p=None, delay=0):
        '''Initialize the based on the passed in option.

        Parameters
        ----------
        value_option (str)
            Different ways to initialize the values. Options:
                'auto', 'all-off'
                    Turn all of the indicators off
                'all-on'
                    Turn all the indicators on
                'random'
                    Randomly assign the indicator with probability `p`
        p : float
            Only required if `value_option` == 'random'
        delay : int
            How many Gibbs steps to delay updating the values
        '''
        # print('in pert ind')
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        self.delay = delay

        for perturbation in self.perturbations:
            perturbation.indicator.set_signal_when_clusters_change(True)
        self.clustering = self.G.perturbations[0].indicator.clustering

        # Set the value
        if not pl.isstr(value_option):
            raise ValueError('`value_option` ({}) must be a str'.format(type(value_option)))
        if value_option in ['all-off', 'auto']:
            value = False
        elif value_option == 'all-on':
            value = True
        elif value_option == 'random':
            if not pl.isfloat(p):
                raise TypeError('`p` ({}) must be a float'.format(type(p)))
            if p < 0 or p > 1:
                raise ValueError('`p` ({}) must be [0,1]'.format(p))
        else:
            raise ValueError('`value_option` ({}) not recognized'.format(value_option))
        for perturbation in self.perturbations:
            for cid in perturbation.clustering.clusters:
                if value_option == 'random':
                    perturbation.indicator.value[cid] = bool(pl.random.bernoulli.sample(p))
                else:
                    perturbation.indicator.value[cid] = value

        # These are for the function `self._make_idx_for_clusters`
        self.ndts_bias = []
        self.n_replicates = self.G.data.n_replicates
        self.n_perturbations = len(self.G.perturbations)
        self.n_dts_for_replicate = self.G.data.n_dts_for_replicate
        self.total_dts = np.sum(self.n_dts_for_replicate)
        self.replicate_bias = np.zeros(self.n_replicates, dtype=int)
        self.n_asvs = len(self.G.data.asvs)
        for ridx in range(1, self.n_replicates):
            self.replicate_bias[ridx] = self.replicate_bias[ridx-1] + \
                self.n_asvs * self.n_dts_for_replicate[ridx - 1]
        for ridx in range(self.G.data.n_replicates):
            self.ndts_bias.append(
                np.arange(0, self.G.data.n_dts_for_replicate[ridx] * self.n_asvs, self.n_asvs))

        s = 'Perturbation indicator initialization results:\n'
        for i, perturbation in enumerate(self.perturbations):
            s += '\tPerturbation {}:\n' \
                '\t\tindicator: {}\n'.format(i, perturbation.indicator.cluster_bool_array())
        logging.info(s)

    def _make_idx_for_clusters(self):
        '''Creates a dictionary that maps the cluster id to the
        rows that correspond to each ASV in the cluster.

        We cannot cast this with numba because it does not support Fortran style
        raveling :(.

        Returns
        -------
        dict: int -> np.ndarray
            Maps the cluster ID to the row indices corresponding to it
        '''
        clusters = [np.asarray(oidxs, dtype=int).reshape(-1,1) \
            for oidxs in self.clustering.tolistoflists()]
        n_dts=self.G.data.n_dts_for_replicate

        d = {}
        cids = self.clustering.order

        for cidx,cid in enumerate(cids):
            a = np.zeros(len(clusters[cidx]) * self.total_dts, dtype=int)
            i = 0
            for ridx in range(self.n_replicates):
                idxs = np.zeros(
                    (len(clusters[cidx]),
                    self.n_dts_for_replicate[ridx]), int)
                idxs = idxs + clusters[cidx]
                idxs = idxs + self.ndts_bias[ridx]
                idxs = idxs + self.replicate_bias[ridx]
                idxs = idxs.ravel('F')
                l = len(idxs)
                a[i:i+l] = idxs
                i += l

            d[cid] = a
        
        if self.G.data.zero_inflation_transition_policy is not None:
            # We need to convert the indices that are meant from no zero inflation to 
            # ones that take into account zero inflation - use the array from 
            # `data.Data._setrows_to_include_zero_inflation`. If the index should be
            # included, then we subtract the number of indexes that are previously off
            # before that index. If it should not be included then we exclude it
            prevoff_arr = self.G.data.off_previously_arr_zero_inflation
            rows_to_include = self.G.data.zero_inflation_transition_policy
            for cid in d:
                arr = d[cid]
                new_arr = np.zeros(len(arr), dtype=int)
                n = 0
                for i, idx in enumerate(arr):
                    if rows_to_include[idx]:
                        new_arr[n] = idx - prevoff_arr[i]
                        n += 1
                new_arr = new_arr[:n]
        return d

    # @profile
    def make_rel_params(self):
        '''We make the parameters needed to update the relative log-likelihod.
        This function is called once at the beginning of the update.

        THIS ASSUMES THAT EACH PERTURBATION CLUSTERS ARE DEFINED BY THE SAME CLUSTERS
            - To make this separate, make a higher level list for each perturbation index
              for each individual perturbation

        Parameters that we create with this function
        --------------------------------------------
        ys : dict (int -> np.ndarray)
            Maps the target cluster id to the observation matrix that it
            corresponds to (only the ASVs in the target cluster). This 
            array already has the growth and self-interactions subtracted
            out:
                $ \frac{log(x_{k+1}) - log(x_{k})}{dt} - a_{1,k} - a_{2,k}x_{k} $
        process_precs : dict (int -> np.ndarray)
            Maps the target cluster id to the vector of the process precision
            that corresponds to the target cluster (only the ASVs in the target
            cluster). This is a 1D array that corresponds to the diagonal of what
            would be the precision matrix.
        interactionXs : dict (int -> np.ndarray)
            Maps the target cluster id to the design matrix for the interactions
            going into that cluster. We pre-index it with the rows and columns
        prior_prec_interaction : dict (int -> np.ndarray)
            Maps the target cluster id to to the diagonal of the prior precision 
            for the interaction values.
        prior_mean_interaction : dict (int -> np.ndarray)
            Maps the target cluster id to to the diagonal of the prior mean 
            for the interaction values.
        prior_ll_ons : np.ndarray
            Prior log likelihood of a positive indicator. These are separate for each
            perturbation.
        prior_ll_offs : np.ndarray
            Prior log likelihood of the negative indicator. These are separate for each
            perturbation.
        priorvar_logdet_diffs : np.ndarray
            This is the prior variance log determinant that we add when the indicator
            is positive. This is different for each perturbation.
        perturbationsXs : dict (int -> np.ndarray)
            Maps the target cluster id to the design matrix that corresponds to 
            the on perturbations of the target clusters. This is preindexed by the 
            rows but not the columns - the columns assume that all of the perturbations
            are on and we index the ones that we want.
        prior_prec_perturbations : np.ndarray
            This is the prior precision of the magnitude for each of the perturbations. Use
            the perturbation index to get the value
        prior_mean_perturbations : np.ndarray
            This is the prior mean of the magnitude for each one of the perturbations. Use
            the perturbation index to get the value
        '''
        row_idxs = self._make_idx_for_clusters()

        # Create ys
        self.ys = {}
        y = self.G.data.construct_lhs(keys=[
            REPRNAMES.SELF_INTERACTION_VALUE, REPRNAMES.GROWTH_VALUE],
            kwargs_dict={REPRNAMES.GROWTH_VALUE:{'with_perturbations': False}})
        for tcid in self.clustering.order:
            self.ys[tcid] = y[row_idxs[tcid], :]

        # Create process_precs
        self.process_precs = {}
        process_prec_diag = self.G[REPRNAMES.PROCESSVAR].prec
        for tcid in self.clustering.order:
            self.process_precs[tcid] = process_prec_diag[row_idxs[tcid]]

        # Make interactionXs
        self.interactionXs = {}
        interactions = self.G[REPRNAMES.INTERACTIONS_OBJ]
        XM_master = self.G.data.design_matrices[REPRNAMES.CLUSTER_INTERACTION_VALUE].toarray()
        for tcid in self.clustering.order:
            cols = []
            for i, interaction in enumerate(interactions.iter_valid()):
                if interaction.target_cid == tcid:
                    if interaction.indicator:
                        cols.append(i)
            cols = np.asarray(cols, dtype=int)
            self.interactionXs[tcid] = pl.util.fast_index(M=XM_master, 
                rows=row_idxs[tcid], cols=cols)

        # Make prior parameters for interactions
        self.prior_prec_interaction = 1/self.G[REPRNAMES.PRIOR_VAR_INTERACTIONS].value
        self.prior_mean_interaction = self.G[REPRNAMES.PRIOR_MEAN_INTERACTIONS].value

        # Make the perturbation parameters
        self.prior_ll_ons = []
        self.prior_ll_offs = []
        self.priorvar_logdet_diffs = []
        self.prior_prec_perturbations = []
        self.prior_mean_perturbations = []

        for perturbation in self.G.perturbations:
            prob_on = perturbation.probability.value
            self.prior_ll_ons.append(np.log(prob_on))
            self.prior_ll_offs.append(np.log(1 - prob_on))
            
            self.priorvar_logdet_diffs.append(
                np.log(perturbation.magnitude.prior.var.value))

            self.prior_prec_perturbations.append( 
                1/perturbation.magnitude.prior.var.value)

            self.prior_mean_perturbations.append(
                perturbation.magnitude.prior.mean.value)

        # Make perturbation matrices
        self.perturbationsXs = {}
        self.G.data.design_matrices[REPRNAMES.PERT_VALUE].M.build(build=True, 
            build_for_neg_ind=True)
        Xpert_master = self.G.data.design_matrices[REPRNAMES.PERT_VALUE].toarray()
        for tcid in self.clustering.order:
            self.perturbationsXs[tcid] = Xpert_master[row_idxs[tcid], :]

        self.n_clusters = len(self.clustering.order)
        self.clustering_order = self.clustering.order

        self.col2pidxcidx = []
        for pidx in range(len(self.perturbations)):
            for cidx in range(len(self.clustering.order)):
                self.col2pidxcidx.append((pidx, cidx))

        self.arr = []
        for perturbation in self.perturbations:
            self.arr = np.append(
                self.arr, 
                perturbation.indicator.cluster_bool_array())
        self.arr = np.asarray(self.arr, dtype=bool)

    # @profile
    def update_relative(self):
        '''Update each perturbation indicator for the given cluster by
        calculating the realtive loglikelihoods of it being on/off as
        supposed to as is. Because this is a relative loglikelihood, we
        only need to take into account the following parameters of the
        model:
            - Only the ASVs in the cluster in question
            - Only the perturbations for that cluster
            - Only the interactions going into the cluster

        Because these matrices are considerably smaller and considered 'dense', we
        do the operations in numpy instead of scipy sparse.

        We permute the order that the indices are updated for more robust mixing
        '''
        if self.sample_iter < self.delay:
            return
        start_time = time.time()

        self.make_rel_params()

        # Iterate over each perturbation indicator variable
        iidxs = npr.permutation(len(self.arr))
        for iidx in iidxs:
            self.update_single_idx_fast(idx=iidx)

        # Set the perturbation indicators from arr
        i = 0
        for perturbation in self.perturbations:
            for cid in self.clustering.order:
                perturbation.indicator.value[cid] = self.arr[i]
                i += 1

        # rebuild the growth design matrix
        self.G.data.design_matrices[REPRNAMES.PERT_VALUE].M.build()
        self.G.data.design_matrices[REPRNAMES.GROWTH_VALUE].build_with_perturbations()
        self._time_taken = time.time() - start_time

    # @profile
    def update_single_idx_fast(self, idx):
        '''Do a Gibbs step for a single cluster
        '''
        pidx, cidx = self.col2pidxcidx[idx]

        prior_ll_on = self.prior_ll_ons[pidx]
        prior_ll_off = self.prior_ll_offs[pidx]

        d_on = self.calculate_relative_marginal_loglikelihood(idx=idx, val=True)
        d_off = self.calculate_relative_marginal_loglikelihood(idx=idx, val=False)

        ll_on = d_on + prior_ll_on
        ll_off = d_off + prior_ll_off
        dd = [ll_off, ll_on]

        # print('\nindicator', idx)
        # print('fast\n\ttotal: {}\n\tbeta_logdet_diff: {}\n\t' \
        #     'priorvar_logdet_diff: {}\n\tbEb_diff: {}\n\t' \
        #     'bEbprior_diff: {}'.format(
        #         ll_on - ll_off,
        #         d_on['beta_logdet'] - d_off['beta_logdet'],
        #         d_on['priorvar_logdet'] - d_off['priorvar_logdet'],
        #         d_on['bEb'] - d_off['bEb'],
        #         d_on['bEbprior'] - d_off['bEbprior']))
        # self.update_single_idx_slow(idx)

        res = bool(sample_categorical_log(dd))
        self.arr[idx] = res

    # @profile
    def calculate_relative_marginal_loglikelihood(self, idx, val):
        '''Calculate the relative marginal loglikelihood of switching the `idx`'th index
        of the perturbation matrix to `val`

        Parameters
        ----------
        idx : int
            This is the index of the indicator we are sampling
        val : bool
            This is the value we are testing it at.

        Returns
        -------
        float
        '''
        # Create and get the data
        self.arr[idx] = val
        pidx, cidx = self.col2pidxcidx[idx]
        tcid = self.clustering_order[cidx]

        y = self.ys[tcid]
        process_prec = self.process_precs[tcid]
        X = self.interactionXs[tcid]

        prior_mean = []
        prior_prec_diag = []
        cols = []
        for temp_pidx in range(len(self.perturbations)):
            col = int(cidx + temp_pidx * self.n_clusters)
            if self.arr[col]:
                cols.append(col)
                prior_mean.append(self.prior_mean_perturbations[temp_pidx])
                prior_prec_diag.append(self.prior_prec_perturbations[temp_pidx])
        Xpert = self.perturbationsXs[tcid][:, cols]

        if Xpert.shape[1] + X.shape[1] == 0:
            # return {
            #     'ret': 0,
            #     'beta_logdet': 0,
            #     'priorvar_logdet': 0,
            #     'bEb': 0,
            #     'bEbprior': 0}
            return 0

        prior_mean = np.append(
            prior_mean,
            np.full(X.shape[1], self.prior_mean_interaction))
        prior_prec_diag = np.append(
            prior_prec_diag,
            np.full(X.shape[1], self.prior_prec_interaction))
        X = np.hstack((Xpert, X))
        prior_prec = np.diag(prior_prec_diag)
        pm = (prior_prec_diag * prior_mean).reshape(-1,1)

        # Do the marginalization
        a = X.T * process_prec
        beta_prec = (a @ X) + prior_prec
        beta_cov = pinv(beta_prec, self)
        beta_mean = beta_cov @ ((a @ y) + pm )

        bEb = (beta_mean.T @ beta_prec @ beta_mean)[0,0]
        try:
            beta_logdet = log_det(beta_cov, self)
        except:
            logging.critical('Crashed in log_det')
            logging.critical('beta_cov:\n{}'.format(beta_cov))
            logging.critical('prior_prec\n{}'.format(prior_prec))
            raise
        
        if val:
            bEbprior = (self.prior_mean_perturbations[pidx]**2) * \
                self.prior_prec_perturbations[pidx]
            priorvar_logdet = self.priorvar_logdet_diffs[pidx]
        else:
            bEbprior = 0
            priorvar_logdet = 0

        ll2 = 0.5 * (beta_logdet - priorvar_logdet)
        ll3 = 0.5 * (bEb - bEbprior)

        # return {
        #     'ret': ll2+ll3,
        #     'beta_logdet': beta_logdet,
        #     'priorvar_logdet': priorvar_logdet,
        #     'bEb': bEb,
        #     'bEbprior': bEbprior}
        return ll2 + ll3        

    # @profile
    def update_slow(self):
        '''Update each cluster indicator variable for the perturbation
        '''
        start_time = time.time()

        if self.sample_iter < self.delay:
            return

        n_clusters = len(self.clustering.order)
        n_perturbations = len(self.perturbations)
        idxs = npr.permutation(int(n_clusters*n_perturbations))
        for idx in idxs:
            self.update_single_idx_slow(idx=idx)

        # rebuild the growth design matrix
        self.G.data.design_matrices[REPRNAMES.PERT_VALUE].M.build()
        self.G.data.design_matrices[REPRNAMES.GROWTH_VALUE].build_with_perturbations()
        self._time_taken = time.time() - start_time

    # @profile
    def update_single_idx_slow(self, idx):
        '''Do a Gibbs step for a single cluster and perturbation

        Parameters
        ----------
        This is the index of the indicator in vectorized form
        '''
        cidx = idx % self.G.data.n_asvs
        cid = self.clustering.order[cidx]
        
        pidx = idx // self.G.data.n_asvs
        perturbation = self.perturbations[pidx]

        d_on = self.calculate_marginal_loglikelihood(cid=cid, val=True,
            perturbation=perturbation)
        d_off = self.calculate_marginal_loglikelihood(cid=cid, val=False,
            perturbation=perturbation)

        prior_ll_on = np.log(perturbation.probability.value)
        prior_ll_off = np.log(1 - perturbation.probability.value)

        ll_on = d_on['ret'] + prior_ll_on
        ll_off = d_off['ret'] + prior_ll_off

        # print('slow\n\ttotal: {}\n\tbeta_logdet_diff: {}\n\t' \
        #     'priorvar_logdet_diff: {}\n\tbEb_diff: {}\n\t' \
        #     'bEbprior_diff: {}'.format(
        #         ll_on - ll_off,
        #         d_on['beta_logdet'] - d_off['beta_logdet'],
        #         d_on['priorvar_logdet'] - d_off['priorvar_logdet'],
        #         d_on['bEb'] - d_off['bEb'],
        #         d_on['bEbprior'] - d_off['bEbprior']))

        dd = [ll_off, ll_on]

        res = bool(sample_categorical_log(dd))
        if perturbation.indicator.value[cid] != res:
            perturbation.indicator.value[cid] = res
            self.G.data.design_matrices[REPRNAMES.PERT_VALUE].build()

    # @profile
    def calculate_marginal_loglikelihood(self, cid, val, perturbation):
        '''Calculate the log marginal likelihood with the perturbations integrated
        out
        '''
        # Set parameters
        perturbation.indicator.value[cid] = val
        self.G.data.design_matrices[REPRNAMES.PERT_VALUE].M.build()

        # Make matrices
        rhs = [REPRNAMES.PERT_VALUE, REPRNAMES.CLUSTER_INTERACTION_VALUE]
        lhs = [REPRNAMES.GROWTH_VALUE, REPRNAMES.SELF_INTERACTION_VALUE]
        X = self.G.data.construct_rhs(keys=rhs)
        y = self.G.data.construct_lhs(keys=lhs, 
            kwargs_dict={REPRNAMES.GROWTH_VALUE:{'with_perturbations': False}})
        
        if X.shape[1] == 0:
            return {
            'ret': 0,
            'beta_logdet': 0,
            'priorvar_logdet': 0,
            'bEb': 0,
            'bEbprior': 0}

        process_prec = self.G[REPRNAMES.PROCESSVAR].build_matrix(cov=False, sparse=False)
        prior_prec = build_prior_covariance(G=self.G, cov=False, order=rhs, sparse=False)
        prior_var = build_prior_covariance(G=self.G, cov=True, order=rhs, sparse=False)
        prior_mean = build_prior_mean(G=self.G, order=rhs, shape=(-1,1))

        # Calculate the posterior
        beta_prec = X.T @ process_prec @ X + prior_prec
        beta_cov = pinv(beta_prec, self)
        beta_mean = beta_cov @ ( X.T @ process_prec @ y + prior_prec @ prior_mean )
        beta_mean = np.asarray(beta_mean).reshape(-1,1)

        # Perform the marginalization
        try:
            beta_logdet = log_det(beta_cov, self)
        except:
            logging.critical('Crashed in log_det')
            logging.critical('beta_cov:\n{}'.format(beta_cov))
            logging.critical('prior_prec\n{}'.format(prior_prec))
            raise
        priorvar_logdet = log_det(prior_var, self)
        ll2 = 0.5 * (beta_logdet - priorvar_logdet)

        a = np.asarray(prior_mean.T @ prior_prec @ prior_mean)[0,0]
        b = np.asarray(beta_mean.T @ beta_prec @ beta_mean)[0,0]
        ll3 = -0.5 * (a  - b)

        return {
            'ret': ll2+ll3,
            'beta_logdet': beta_logdet,
            'priorvar_logdet': priorvar_logdet,
            'bEb': b,
            'bEbprior': a}

    def total_on(self):
        n = 0
        for perturbation in self.perturbations:
            n += perturbation.indicator.num_on_clusters()
        return n


class PriorVarPerturbations(pl.Variable):
    '''Agglomerates the prior variances of the magnitudes for the perturbations.

    All perturbations get the same hyperparameters
    '''
    def __init__(self, **kwargs):
        kwargs['name'] = STRNAMES.PRIOR_VAR_PERT
        pl.Variable.__init__(self, **kwargs)
        self.perturbations = self.G.perturbations

        if self.perturbations is None:
            raise TypeError('Only instantiate this object if there are perturbations')

    def __str__(self):
        s = 'Perturbation Magnitude Prior Variances'
        for perturbation in self.perturbations:
            s += '\n\tperturbation {}: {}'.format(
                perturbation.name,
                perturbation.magnitude.prior.var.value)
        return s

    @property
    def sample_iter(self):
        return self.perturbations[0].magnitude.prior.var.sample_iter

    def initialize(self, **kwargs):
        '''Every prior variance on the perturbations gets the same hyperparameters
        '''
        for perturbation in self.perturbations:
            perturbation.magnitude.prior.var.initialize(**kwargs)

    def update(self):
        for perturbation in self.perturbations:
            perturbation.magnitude.prior.var.update()

    def set_trace(self, *args, **kwargs):
        for perturbation in self.perturbations:
            perturbation.magnitude.prior.var.set_trace(*args, **kwargs)

    def add_trace(self, *args, **kwargs):
        for perturbation in self.perturbations:
            perturbation.magnitude.prior.var.add_trace(*args, **kwargs)

    def get_single_value_of_perts(self):
        '''Get the variance for each perturbation
        '''
        return np.asarray([p.magnitude.prior.var.value for p in self.perturbations])

    def diag(self, only_pos_ind=True):
        '''Return the diagonal of the prior variances stacked up in order

        Parameters
        ----------
        only_pos_ind : bool
            If True, only put in the values for the positively indicated clusters
            for each perturbation

        Returns
        -------
        np.ndarray
        '''
        ret = []
        for perturbation in self.perturbations:
            if only_pos_ind:
                n = perturbation.indicator.num_on_clusters()
            else:
                n = len(perturbation.clustering)
            ret = np.append(
                ret,
                np.ones(n, dtype=float)*perturbation.magnitude.prior.var.value)
        return ret


class PriorVarPerturbationSingle(pl.variables.SICS):
    '''This is the posterior of the prior variance of regression coefficients
    for the interaction (off diagonal) variables
    '''
    def __init__(self, prior, perturbation, value=None, **kwargs):

        kwargs['name'] = STRNAMES.PRIOR_VAR_PERT + '_' + perturbation.name
        pl.variables.SICS.__init__(self, value=value, dtype=float, **kwargs)
        self.add_prior(prior)
        self.perturbation = perturbation

    def initialize(self, value_option, dof_option, scale_option, value=None,
        dof=None, scale=None, delay=0):
        '''Initialize the hyperparameters of the perturbation prior variance based on the
        passed in option

        Parameters
        ----------
        value_option : str
            - Initialize the value based on the specified option
            - Options
                'manual'
                    Set the value manually, `value` must also be specified
                'auto', 'prior-mean'
                    Set the value to the mean of the prior
                'tight'
                    value = 10^2
                'diffuse'
                    value = 10^4
        scale_option : str
            Initialize the scale of the prior
            Options
                'manual'
                    Set the value manually, `scale` must also be specified
                'auto', 'diffuse'
                    Set so that the mean of the distribution is 10^4
                'tight'
                    Set so that the mean of the distribution is 10^2
        dof_option : str
            Initialize the dof of the parameter
            Options:
                'manual': Set the value with the parameter `dof`
                'diffuse': Set the value to 2.5
                'strong': Set the value to the expected number of interactions
                'auto': Set to diffuse
        dof, scale : int, float
            User specified values
            Only necessary if  any of the options are 'manual'
        '''
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        self.delay = delay

        if not pl.isstr(dof_option):
            raise TypeError('`dof_option` ({}) must be a str'.format(type(dof_option)))
        if dof_option == 'manual':
            if not pl.isnumeric(dof):
                raise TypeError('`dof` ({}) must be a numeric'.format(type(dof)))
            if dof < 0:
                raise ValueError('`dof` ({}) must be > 0 for it to be a valid prior'.format(shape))
        elif dof_option in ['diffuse', 'auto']:
            dof = 2.5
        elif dof_option == 'strong':
            dof = expected_n_clusters(G=self.G)
        else:
            raise ValueError('`dof_option` ({}) not recognized'.format(dof_option))
        self.prior.dof.override_value(dof)

        if not pl.isstr(scale_option):
            raise TypeError('`scale_option` ({}) must be a str'.format(type(scale_option)))
        if scale_option == 'manual':
            if not pl.isnumeric(scale):
                raise TypeError('`scale` ({}) must be a numeric'.format(type(scale)))
            if scale < 0:
                raise ValueError('`scale` ({}) must be > 0 for it to be a valid prior'.format(scale))
        elif scale_option in ['auto', 'diffuse']:
            # Calculate the mean to be 10
            scale = 1e4 * (self.prior.dof.value - 2) / self.prior.dof.value
        elif scale_option == 'tight':
            scale = 100 * (self.prior.dof.value - 2) / self.prior.dof.value
        else:
            raise ValueError('`scale_option` ({}) not recognized'.format(scale_option))
        self.prior.scale.override_value(scale)

        if not pl.isstr(value_option):
            raise TypeError('`value_option` ({}) must be a str'.format(type(value_option)))
        if value_option == 'manual':
            if not pl.isnumeric(value):
                raise ValueError('`value` ({}) must be numeric (float,int)'.format(value.__class__))
            self.value = value
        elif value_option in ['auto', 'prior-mean']:
            self.value = self.prior.mean()
        elif value_option == 'diffuse':
            self.value = 1e4
        elif value_option == 'tight':
            self.value = 1e2
        else:
            raise ValueError('`value_option` ({}) not recognized'.format(value_option))

        logging.info('Prior Variance Interactions initialization results:\n' \
            '\tprior dof: {}\n' \
            '\tprior scale: {}\n' \
            '\tvalue: {}'.format(
                self.prior.dof.value, self.prior.scale.value, self.value))

    # @profile
    def update(self):
        '''Calculate the posterior of the prior variance
        '''
        if self.sample_iter < self.delay:
            return

        x = self.perturbation.cluster_array(only_pos_ind=True)
        mu = self.perturbation.magnitude.prior.mean.value

        se = np.sum(np.square(x - mu))
        n = len(x)

        self.dof.value = self.prior.dof.value + n
        self.scale.value = ((self.prior.scale.value * self.prior.dof.value) + \
           se)/self.dof.value
        self.sample()


class PriorMeanPerturbations(pl.Variable):
    '''Agglomerates the prior variances of the magnitudes for the perturbations.

    All perturbations get the same hyperparameters
    '''
    def __init__(self, **kwargs):
        kwargs['name'] = STRNAMES.PRIOR_MEAN_PERT
        pl.Variable.__init__(self, **kwargs)
        self.perturbations = self.G.perturbations

        if self.perturbations is None:
            raise TypeError('Only instantiate this object if there are perturbations')

    def __str__(self):
        s = 'Perturbation Magnitude Prior Means'
        for perturbation in self.perturbations:
            s += '\n\tperturbation {}: {}'.format(
                perturbation.name,
                perturbation.magnitude.prior.mean.value)
        return s

    @property
    def sample_iter(self):
        return self.perturbations[0].magnitude.prior.mean.sample_iter

    def initialize(self, **kwargs):
        '''Every prior variance on the perturbations gets the same hyperparameters
        '''
        for perturbation in self.perturbations:
            perturbation.magnitude.prior.mean.initialize(**kwargs)

    def update(self):
        for perturbation in self.perturbations:
            perturbation.magnitude.prior.mean.update()

    def set_trace(self, *args, **kwargs):
        for perturbation in self.perturbations:
            perturbation.magnitude.prior.mean.set_trace(*args, **kwargs)

    def add_trace(self, *args, **kwargs):
        for perturbation in self.perturbations:
            perturbation.magnitude.prior.mean.add_trace(*args, **kwargs)

    def get_single_value_of_perts(self):
        '''Get the variance for each perturbation
        '''
        return np.asarray([p.magnitude.prior.mean.value for p in self.perturbations])

    def toarray(self, only_pos_ind=True):
        '''Return the diagonal of the prior variances stacked up in order

        Parameters
        ----------
        only_pos_ind : bool
            If True, only put in the values for the positively indicated clusters
            for each perturbation

        Returns
        -------
        np.ndarray
        '''
        ret = []
        for perturbation in self.perturbations:
            if only_pos_ind:
                n = perturbation.indicator.num_on_clusters()
            else:
                n = len(perturbation.clustering)
            ret = np.append(
                ret,
                np.ones(n, dtype=float)*perturbation.magnitude.prior.var.value)
        return ret


class PriorMeanPerturbationSingle(pl.variables.Normal):
    
    def __init__(self, prior, perturbation, **kwargs):

        kwargs['name'] = STRNAMES.PRIOR_MEAN_PERT + '_' + perturbation.name
        pl.variables.Normal.__init__(self, mean=None, var=None, dtype=float, **kwargs)
        self.add_prior(prior)
        self.perturbation = perturbation

    def initialize(self, value_option, mean_option, var_option, value=None,
        mean=None, var=None, delay=0):
        '''Initialize the hyperparameters

        Parameters
        ----------
        value_option : str
            How to set the value. Options:
                'zero'
                    Set to zero
                'prior-mean', 'auto'
                    Set to the mean of the prior
                'manual'
                    Specify with the `value` parameter
        mean_option : str
            How to set the mean of the prior
                'zero', 'auto'
                    Set to zero
                'manual'
                    Set with the `mean` parameter
        var_option : str
            'diffuse', 'auto'
                Variance is set to 10e4
            'tight'
                Variance is set to 1e2
            'manual'
                Set with the `var` parameter
        value, mean, var : float
            These are only necessary if we specify manual for any of the other 
            options
        delay : int
            How much to delay the start of the update during inference
        '''
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        self.delay = delay

        # Set the mean
        if not pl.isstr(mean_option):
            raise TypeError('`mean_option` ({}) must be a str'.format(type(mean_option)))
        if mean_option == 'manual':
            if not pl.isnumeric(mean):
                raise TypeError('`mean` ({}) must be a numeric'.format(type(mean)))
        elif mean_option in ['zero', 'auto']:
            mean = 0
        else:
            raise ValueError('`mean_option` ({}) not recognized'.format(mean_option))
        self.prior.mean.override_value(mean)

        # Set the variance
        if not pl.isstr(var_option):
            raise TypeError('`var_option` ({}) must be a str'.format(type(var_option)))
        if var_option == 'manual':
            if not pl.isnumeric(var):
                raise TypeError('`var` ({}) must be a numeric'.format(type(var)))
            if var <= 0:
                raise ValueError('`var` ({}) must be positive'.format(var))
        elif var_option in ['diffuse', 'auto']:
            var = 1e4
        elif var_option == 'tight':
            var = 1e2
        else:
            raise ValueError('`var_option` ({}) not recognized'.format(var_option))
        self.prior.var.override_value(var)

        # Set the value
        if not pl.isstr(value_option):
            raise TypeError('`value_option` ({}) must be a str'.format(type(value_option)))
        if value_option == 'manual':
            if not pl.isnumeric(value):
                raise TypeError('`value` ({}) must be a numeric'.format(type(value)))
        elif value_option in ['prior-mean', 'auto']:
            value = self.prior.mean.value
        elif value_option == 'zero':
            value = 0
        else:
            raise ValueError('`value_option` ({}) not recognized'.format(value_option))
        self.value = value
        
    def update(self):
        '''Update using a Gibbs update
        '''
        if self.sample_iter < self.delay:
            return

        x = self.perturbation.cluster_array(only_pos_ind=True)
        prec = 1/self.perturbation.magnitude.prior.var.value

        prior_prec = 1/self.prior.var.value
        prior_mean = self.prior.mean.value

        self.var.value = 1/(prior_prec + (len(x)*prec))
        self.mean.value = self.var.value * ((prior_mean * prior_prec) + (np.sum(x)*prec))
        self.sample()


# qPCR variance (NOTE THIS IS NOT USED IN THE MODEL)
# --------------------------------------------------
class _qPCRBase(pl.Variable):
    '''Base class for qPCR measurements
    '''
    def __init__(self, L, **kwargs):

        pl.Variable.__init__(self, **kwargs)
        self._sample_iter = 0
        self.L = L
        self.n_replicates = self.G.data.n_replicates
        self.value = []

    def __getitem__(self, key):
        return self.value[key]

    def __str__(self):
        # Make them into an array?
        try:
            s = ''
            for a in self.value:
                s += str(a) + '\n'
        except:
            s = 'not set'
        return s

    def update(self):
        '''Update each of the qPCR variances
        '''
        for a in self.value:
            a.update()
        self._sample_iter += 1

    def initialize(self, **kwargs):
        '''Every variance gets the same initialization
        '''
        for a in self.value:
            a.initialize(**kwargs)

    def add_trace(self, *args, **kwargs):
        for a in self.value:
            a.add_trace(*args, **kwargs)

    def set_trace(self, *args, **kwargs):
        for a in self.value:
            a.set_trace(*args, **kwargs)


class _qPCRPriorAggVar(_qPCRBase):
    '''Base class for `qPCRDegsOfFreedoms` and `qPCRScales`
    '''
    def __init__(self, L, child,**kwargs):
        _qPCRBase.__init__(self, L, **kwargs)
        for l in range(self.L):
            self.value.append(
                child(L=L, l=l, **kwargs))

    def add_qpcr_measurement(self, ridx, tidx, l):
        '''Add a qPCR measurement for subject `ridx` at time index
        `tidx` to qPCR set `l`

        Parameters
        ----------
        ridx : int
            Subject index
        tidx : int
            Time index
        l : int
            qPCR set index
        '''
        if not pl.isint(ridx):
            raise TypeError('`ridx` ({}) must be an int'.format(type(ridx)))
        if ridx >= self.G.data.n_replicates:
            raise ValueError('`ridx` ({}) out of range ({})'.format(ridx, 
                self.G.data.n_replicates))
        if not pl.isint(tidx):
            raise TypeError('`tidx` ({}) must be an int'.format(type(tidx)))
        if tidx >= len(self.G.data.given_timepoints[ridx]):
            raise ValueError('`tidx` ({}) out of range ({})'.format(tidx, 
                len(self.G.data.given_timepoints[ridx])))
        if not pl.isint(l):
            raise TypeError('`l` ({}) must be an int'.format(type(l)))
        if l >= self.L:
            raise ValueError('`l` ({}) out of range ({})'.format(tidx, 
                self.L))
        self.value[l].add_qpcr_measurement(ridx=ridx, tidx=tidx)

    def set_shape(self):
        for a in self.value:
            a.set_shape()


class qPCRVariances(_qPCRBase):
    '''Aggregation class for qPCR variance for a set of qPCR variances. 
    The qPCR variances are 

    Parameters
    ----------
    L : int
        How many qPCR variance groupings there are
    '''
    def __init__(self, **kwargs):
        kwargs['name'] = STRNAMES.QPCR_VARIANCES
        _qPCRBase.__init__(self, **kwargs)
        self.G.data.qpcr_variances = self
        
        for ridx in range(self.n_replicates):
            self.value.append( 
                qPCRVarianceReplicate(ridx=ridx, **kwargs))

    def add_qpcr_measurement(self, ridx, tidx, l):
        '''Add a qPCR measurement for subject `ridx` at time index
        `tidx` to qPCR set `l`

        Parameters
        ----------
        ridx : int
            Subject index
        tidx : int
            Time index
        l : int
            qPCR set index
        '''
        if not pl.isint(ridx):
            raise TypeError('`ridx` ({}) must be an int'.format(type(ridx)))
        if ridx >= self.G.data.n_replicates:
            raise ValueError('`ridx` ({}) out of range ({})'.format(ridx, 
                self.G.data.n_replicates))
        if not pl.isint(tidx):
            raise TypeError('`tidx` ({}) must be an int'.format(type(tidx)))
        if tidx >= len(self.G.data.given_timepoints[ridx]):
            raise ValueError('`tidx` ({}) out of range ({})'.format(tidx, 
                len(self.G.data.given_timepoints[ridx])))
        if not pl.isint(l):
            raise TypeError('`l` ({}) must be an int'.format(type(l)))
        if l >= self.L:
            raise ValueError('`l` ({}) out of range ({})'.format(tidx, 
                self.L))
        self.value[ridx].add_qpcr_measurement(tidx=tidx, l=l)        


class qPCRVarianceReplicate(pl.variables.SICS):
    '''Posterior for a set of single qPCR variances for replicate `ridx`

    Parameters
    ----------
    ridx : int
        Which subject replicate index this set of qPCR variances belongs
        to.
    L : int
        How many qPCR variance groupings there are
    '''
    def __init__(self, ridx, L, **kwargs):
        self.ridx = ridx
        self.L = L
        kwargs['name'] = STRNAMES.QPCR_VARIANCES + '_{}'.format(ridx)
        pl.variables.SICS.__init__(self, **kwargs)
        self.priors_idx = np.full(len(self.G.data.given_timepoints[ridx]), -1, dtype=int)
        self.set_value_shape(shape=(len(self.G.data.given_timepoints[ridx]),))

    def initialize(self, value_option, value=None, inflated=None):
        '''Initialize the values. We do not set any hyperparameters because those are
        set in their own classes.

        Parameters
        ----------
        value_option : str
            How to initialize the variances
                'empirical', 'auto'
                    Set to the empirical variance of the respective measurements
                'inflated'
                    Set to an inflated value of the empirical variance.
                'manual'
                    Set the values manually
        value : float, np.ndarray(float)
            If float, set all the values to the same number. If array then set the 
            values to each of the parameters
        inflated : float, None
            Necessary if `value_option` == 'inflated'
        '''
        if not pl.isstr(value_option):
            raise TypeError('`value_option` ({}) must be a str'.format(type(value_option)))
        if value_option in ['empirical', 'auto']:
            self.value = np.zeros(len(self.G.data.qpcr[self.ridx]), dtype=float)
            for idx, t in enumerate(self.G.data.qpcr[self.ridx]):
                self.value[idx] = np.var(self.G.data.qpcr[self.ridx][t].log_data)

        elif value_option == 'inflated':
            if not pl.isnumeric(inflated):
                raise TypeError('`inflated` ({}) must be a numeric'.format(type(inflated)))
            if inflated < 0:
                raise ValueError('`inflated` ({}) must be positive'.format(inflated))
            # Set each variance by the empirical variance * inflated
            self.value = np.zeros(len(self.G.data.qpcr[self.ridx]), dtype=float)
            for idx, t in enumerate(self.G.data.qpcr[self.ridx]):
                self.value[idx] = np.var(self.G.data.qpcr[self.ridx][t].log_data) * inflated

        elif value_option == 'manual':
            raise NotImplementedError('Need to implement')
        else:
            raise ValueError('`value_option` ({}) not recognized'.format(value_option))

        # Set the qPCR measurements
        self.qpcr_measurements = []
        for tidx, t in enumerate(self.G.data.given_timepoints[self.ridx]):
            self.qpcr_measurements.append(self.G.data.qpcr[self.ridx][t].log_data)

    def update(self):

        prior_dofs = []
        prior_scales = []
        for l in range(self.L):
            prior_dofs.append(self.G[REPRNAMES.QPCR_DOFS].value[l].value)
            prior_scales.append(self.G[REPRNAMES.QPCR_SCALES].value[l].value)

        for tidx in range(len(self.priors_idx)):
            t = self.G.data.given_timepoints[self.ridx][tidx]
            l = self.priors_idx[tidx]
            prior_dof = prior_dofs[l]
            prior_scale = prior_scales[l]

            # qPCR measurements (these are already in log space)
            values = self.qpcr_measurements[tidx]

            # Current mean is the log of the sum of latent abundance
            tidx_in_arr = self.G.data.timepoint2index[self.ridx][t]
            mean = np.log(np.sum(self.G.data.data[self.ridx][:, tidx_in_arr]))

            # Calculate the residual sum
            resid_sum = np.sum(np.square(values - mean))

            # posterior
            dof = prior_dof + len(values)
            scale = ((prior_scale * prior_dof) + resid_sum)/dof
            self.value[tidx] = pl.random.sics.sample(dof, scale)

    def add_qpcr_measurement(self, tidx, l):
        '''Add qPCR measurement for subject index `ridx` and time index `tidx`

        Parameters
        ----------
        ridx : int
            Subject index
        tidx : int
            Time index
        '''
        if not pl.isint(tidx):
            raise TypeError('`tidx` ({}) must be an int'.format(type(tidx)))
        if tidx >= len(self.G.data.given_timepoints[self.ridx]):
            raise ValueError('`tidx` ({}) out of range ({})'.format(tidx, 
                len(self.G.data.given_timepoints[self.ridx])))
        if not pl.isint(l):
            raise TypeError('`l` ({}) must be an int'.format(type(l)))
        self.priors_idx[tidx] = l


class qPCRDegsOfFreedoms(_qPCRPriorAggVar):
    '''Aggregation class for a degree of freedom parameter of qPCR variance

    Parameters
    ----------
    L : int
        How many qPCR variance groupings there are
    '''
    def __init__(self, L, **kwargs):
        kwargs['name'] = STRNAMES.QPCR_DOFS
        _qPCRPriorAggVar.__init__(self, L=L, child=qPCRDegsOfFreedomL, **kwargs)


class qPCRDegsOfFreedomL(pl.variables.Uniform):
    '''Posterior for a single qPCR degrees of freedom parameter for a SICS set
    
    Parameters
    ----------
    L : int
        How many qPCR variance groupings there are
    l : int
        Which specific grouping this hyperprior is
    '''
    def __init__(self, L, l, **kwargs):

        self.L = L
        self.l = l
        kwargs['name'] = STRNAMES.QPCR_DOFS + '_{}'.format(l)
        pl.variables.Uniform.__init__(self, **kwargs)

        self.data_locs = []
        self.proposal = pl.variables.TruncatedNormal(mean=None, var=None, value=None)

    def __str__(self):
        # If this fails, it is because we are dividing by 0 sampler_iter
        # If which case we just return the value 
        try:
            s = 'Value: {}, Acceptance rate: {}'.format(
                self.value, np.mean(self.acceptances[
                    np.max([self.sample_iter-50, 0]):self.sample_iter]))
        except:
            s = str(self.value)
        return s

    def set_shape(self):
        '''Set the shape of the array (how many qPCR variances this is a prior for)
        '''
        self.set_value_shape(shape=(len(self.data_locs), ))

    def add_qpcr_measurement(self, ridx, tidx):
        '''Add the qPCR measurement for subject index `ridx` and time index
        `tidx` to 
        '''
        self.data_locs.append((ridx, tidx))

    def initialize(self, value_option, low_option, high_option, proposal_option, 
        target_acceptance_rate, tune, end_tune, value=None, low=None, high=None, 
        proposal_var=None, delay=0):
        '''Initialize the values and hyperparameters. The proposal truncation is 
        always set to the same as the parameterization of the prior.

        Parameters
        ----------
        value_option : str
            How to initialize the value. Options:
                'auto', 'diffuse'
                    Set the value to 2.5
                'strong'
                    Set to be 50% of the data
                'manual'
                    `value` must also be specified
        low_option : str
            How to set the low parameter of the prior
            'auto', 'valid'
                Set to 2 so that the prior stays proper during inference
            'zero'
                Set to 0
            'manual'
                Specify the value with the parameter `low`
        high_option : str
            How to set the high parameter of the prior
            'auto', 'med'
                Set to 10 X the maximum in the set
            'high'
                Set to 100 X the maximum in the set
            'low'
                Set to 1 X the maximum in the set
            'manual'
                Set the value with the parameter `high`
        proposal_option : str
            How to initialize the proposal variance:
                'auto'
                    mean**2 / 100
                'manual'
                    `proposal_var` must also be supplied
        '''
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        self.delay = delay

        self.qpcr_data = []
        for ridx, tidx in self.data_locs:
            t = self.G.data.given_timepoints[ridx][tidx]
            self.qpcr_data = np.append(self.qpcr_data, 
                self.G.data.qpcr[ridx][t].log_data)
        
        # Set the prior low
        if not pl.isstr(low_option):
            raise TypeError('`low_option` ({}) must be a str'.format(type(low_option)))
        if low_option == 'manual':
            if not pl.isnumeric(low):
                raise TypeError('`low` ({}) must be a numeric'.format(type(low)))
            if low < 2:
                raise ValueError('`low` ({}) must be >= 2'.format(low))
        elif low_option in ['valid', 'auto']:
            low = 2
        elif low_option == 'zero':
            low = 0
        else:
            raise ValueError('`low_option` ({}) not recognized'.format(low_option))
        self.prior.low.override_value(low)  

        # Set the prior high
        if not pl.isstr(high_option):
            raise TypeError('`high_option` ({}) must be a str'.format(type(high_option)))
        if high_option == 'manual':
            if not pl.isnumeric(high):
                raise TypeError('`high` ({}) must be a numeric'.format(type(high)))
        elif high_option in ['med', 'auto']:
            high = 10 * len(self.data_locs)
        elif high_option == 'low':
            high = len(self.data_locs)
        elif high_option == 'high':
            high = 100 * len(self.data_locs)
        else:
            raise ValueError('`high_option` ({}) not recognized'.format(high_option))
        if high < self.prior.low.value:
            raise ValueError('`high` ({}) must be >= low ({})'.format(high, 
                self.prior.low.value))
        self.prior.high.override_value(high)

        # Set the value
        if not pl.isstr(value_option):
            raise TypeError('`value_option` ({}) must be a str'.format(type(value_option)))
        if value_option == 'manual':
            if not pl.isnumeric(value):
                raise TypeError('`value` ({}) must be a numeric'.format(type(value)))
        elif value_option in ['auto', 'diffuse']:
            value = 2.5
        elif value_option == 'strong':
            value = len(self.data_locs)
        else:
            raise ValueError('`value_option` ({}) not recognized'.format(value_option))
        self.value = value
        if self.value <= self.prior.low.value or self.value >= self.prior.high.value:
            raise ValueError('`value` ({}) out of range ({})'.format(self.value))

        # Set the propsal parameters
        if pl.isstr(target_acceptance_rate):
            if target_acceptance_rate in ['optimal', 'auto']:
                target_acceptance_rate = 0.44
            else:
                raise ValueError('`target_acceptance_rate` ({}) not recognized'.format(
                    target_acceptance_rate))
        elif pl.isfloat(target_acceptance_rate):
            if target_acceptance_rate < 0 or target_acceptance_rate > 1:
                raise ValueError('`target_acceptance_rate` ({}) out of range'.format(
                    target_acceptance_rate))
        else:
            raise TypeError('`target_acceptance_rate` ({}) type not recognized'.format(
                type(target_acceptance_rate)))
        self.target_acceptance_rate = target_acceptance_rate
        
        if pl.isstr(tune):
            if tune in ['auto']:
                tune = 50
            else:
                raise ValueError('`tune` ({}) not recognized'.format(tune))
        elif pl.isint(tune):
            if tune < 0:
                raise ValueError('`tune` ({}) must be > 0'.format(
                    tune))
        else:
            raise TypeError('`tune` ({}) type not recognized'.format(type(tune)))
        self.tune = tune

        if pl.isstr(end_tune):
            if end_tune in ['auto', 'half-burnin']:
                end_tune = int(self.G.inference.burnin/2)
            else:
                raise ValueError('`tune` ({}) not recognized'.format(end_tune))
        elif pl.isint(end_tune):
            if end_tune < 0 or end_tune > self.G.inference.burnin:
                raise ValueError('`end_tune` ({}) out of range (0, {})'.format(
                    end_tune, self.G.inference.burnin))
        else:
            raise TypeError('`end_tune` ({}) type not recognized'.format(type(end_tune)))
        self.end_tune = end_tune

        # Set the proposal variance
        if not pl.isstr(proposal_option):
            raise TypeError('`proposal_option` ({}) must be a str'.format(
                type(proposal_option)))
        elif proposal_option == 'manual':
            if not pl.isnumeric(proposal_var):
                raise TypeError('`proposal_var` ({}) must be a numeric'.format(
                    type(proposal_var)))
            if proposal_var <= 0:
                raise ValueError('`proposal_var` ({}) not proper'.format(proposal_var))
        elif proposal_option in ['auto']:
            proposal_var = (self.value ** 2)/10
        else:
            raise ValueError('`proposal_option` ({}) not recognized'.format(
                proposal_option))
        self.proposal.var.value = proposal_var
        self.proposal.low = self.prior.low.value
        self.proposal.high = self.prior.high.value

    def update_var(self):
        '''Update the variance of the proposal
        '''
        if self.sample_iter == 0:
            self.temp_acceptances = 0
            self.acceptances = np.zeros(self.G.inference.n_samples, dtype=bool)
        
        elif self.sample_iter > self.end_tune:
            # Don't do any more updates
            return
        
        elif self.sample_iter % self.tune == 0:
            # Update var
            acceptance_rate = self.temp_acceptances / self.tune
            if acceptance_rate > self.target_acceptance_rate:
                self.proposal.var.value *= 1.5
            else:
                self.proposal.var.value /= 1.5
            self.temp_acceptances = 0

    def update(self):
        '''First we update the proposal (if necessary) and then we do a MH step
        '''
        self.update_var()
        proposal_std = np.sqrt(self.proposal.var.value)

        # Get the data
        xs = []
        for ridx, tidx in self.data_locs:
            xs.append(self.G[REPRNAMES.QPCR_VARIANCES].value[ridx].value[tidx])

        # Get the scale
        scale = self.G[REPRNAMES.QPCR_SCALES].value[self.l].value

        # Propose a new value for the dof
        prev_dof = self.value
        self.proposal.mean.value = self.value
        new_dof = self.proposal.sample()

        if new_dof < self.prior.low.value or new_dof > self.prior.high.value:
            # Automatic reject
            self.value = prev_dof
            return

        # Calculate the target distribution log likelihood
        prev_target_ll = 0
        for x in xs:
            prev_target_ll += pl.random.sics.logpdf(value=x,
                scale=scale, dof=prev_dof)
        new_target_ll = 0
        for x in xs:
            new_target_ll += pl.random.sics.logpdf(value=x,
                scale=scale, dof=new_dof)

        # Normalize by the loglikelihood of the proposal
        prev_prop_ll = pl.random.truncnormal.logpdf(
            value=prev_dof, mean=new_dof, std=proposal_std,
            low=self.proposal.low, high=self.proposal.high)
        new_prop_ll = pl.random.truncnormal.logpdf(
            value=new_dof, mean=prev_dof, std=proposal_std,
            low=self.proposal.low, high=self.proposal.high)

        # Accept or reject
        r = (new_target_ll - prev_prop_ll) - \
            (prev_target_ll - new_prop_ll)
        u = np.log(pl.random.misc.fast_sample_standard_uniform())

        # print('\n\n\n{} prior_mean\n----------'.format(self.child_name))
        # print('x', x)
        # print('prev_dof', prev_dof)
        # print('prev_target_ll', prev_target_ll)
        # print('prev_prop_ll', prev_prop_ll)
        # print('new dof', new_dof)
        # print('new_target_ll', new_target_ll)
        # print('new_prop_ll', new_prop_ll)
        # print('\nr', r, u)
            
        if r >= u:
            self.acceptances[self.sample_iter] = True
            self.value = new_dof
            self.temp_acceptances += 1
        else:
            self.value = prev_dof


class qPCRScales(_qPCRPriorAggVar):
    '''Aggregation class for a scale parameter of qPCR variance

    Parameters
    ----------
    L : int
        How many qPCR variance groupings there are
    '''
    def __init__(self, L, **kwargs):
        kwargs['name'] = STRNAMES.QPCR_SCALES
        _qPCRPriorAggVar.__init__(self, L=L, child=qPCRScaleL, **kwargs)


class qPCRScaleL(pl.variables.SICS):
    '''Posterior for a single qPCR scale set
    '''
    def __init__(self, L, l, **kwargs):

        self.L = L
        self.l = l
        kwargs['name'] = STRNAMES.QPCR_SCALES + '_{}'.format(l)
        pl.variables.Uniform.__init__(self, **kwargs)

        self.data_locs = []
        self.proposal = pl.variables.TruncatedNormal(mean=None, var=None, value=None)

    def __str__(self):
        # If this fails, it is because we are dividing by 0 sampler_iter
        # If which case we just return the value 
        try:
            s = 'Value: {}, Acceptance rate: {}'.format(
                self.value, np.mean(self.acceptances[
                    np.max([self.sample_iter-50, 0]):self.sample_iter]))
        except:
            s = str(self.value)
        return s

    def set_shape(self):
        '''Set the shape of the array (how many qPCR variances this is a prior for)
        '''
        self.set_value_shape(shape=(len(self.data_locs), ))

    def add_qpcr_measurement(self, ridx, tidx):
        '''Add the qPCR measurement for subject index `ridx` and time index
        `tidx` to 
        '''
        self.data_locs.append((ridx, tidx))

    def initialize(self, value_option, scale_option, dof_option, proposal_option, 
        target_acceptance_rate, tune, end_tune, value=None, dof=None, scale=None,
        proposal_var=None, delay=0):
        '''Initialize the values and hyperparameters

        Parameters
        ----------
        value_option : str
            How to initialize the value. Options:
                'auto', 'prior-mean'
                    Set to the prior mean
                'manual'
                    `value` must also be specified
        dof_option : str
            How to set the prior dof
                'auto', 'diffuse'
                    2.5
                'manual':
                    set with `dof`
        scale_option : str
            How to set the prior scale
                'empirical', 'auto'
                    Set to the variance of the data assigned to the set
                'manual'
                    Set with the parameter `scale`
        proposal_option : str
            How to initialize the proposal variance:
                'auto'
                    mean**2 / 100
                'manual'
                    `proposal_var` must also be supplied     
        '''
        if not pl.isint(delay):
            raise TypeError('`delay` ({}) must be an int'.format(type(delay)))
        if delay < 0:
            raise ValueError('`delay` ({}) must be >= 0'.format(delay))
        self.delay = delay

        self.qpcr_data = []
        for ridx, tidx in self.data_locs:
            t = self.G.data.given_timepoints[ridx][tidx]
            self.qpcr_data = np.append(self.qpcr_data, 
                self.G.data.qpcr[ridx][t].log_data)

        # Set the prior dof
        if not pl.isstr(dof_option):
            raise TypeError('`dof_option` ({}) must be a str'.format(type(dof_option)))
        if dof_option == 'manual':
            if not pl.isnumeric(dof):
                raise TypeError('`dof` ({}) must be a numeric'.format(type(dof)))
            if dof < 2:
                raise ValueError('`dof` ({}) must be >= 2'.format(dof))
        elif dof_option in ['diffuse', 'auto']:
            dof = 2.5
        else:
            raise ValueError('`dof_option` ({}) not recognized'.format(dof_option))
        if dof < 2:
            raise ValueError('`dof` ({}) must be strictly larger than 2 to be a proper' \
                ' prior'.format(dof))
        self.prior.dof.override_value(dof)

        # Set the prior scale
        if not pl.isstr(scale_option):
            raise TypeError('`scale_option` ({}) must be a str'.format(type(scale_option)))
        if scale_option == 'manual':
            if not pl.isnumeric(scale):
                raise TypeError('`scale` ({}) must be a numeric'.format(type(scale)))
            if scale <= 0:
                raise ValueError('`scale` ({}) must be positive'.format(scale))
        elif scale_option in ['auto', 'empirical']:
            v = np.var(self.qpcr_data)
            scale = v * (self.prior.dof.value - 2) / self.prior.dof.value
        else:
            raise ValueError('`scale_option` ({}) not recognized'.format(scale_option))
        self.prior.scale.override_value(scale)

        # Set the value
        if not pl.isstr(value_option):
            raise TypeError('`value_option` ({}) must be a str'.format(type(value_option)))
        if value_option == 'manual':
            if not pl.isnumeric(value):
                raise TypeError('`value` ({}) must be a numeric'.format(type(value)))
        elif value_option in ['auto', 'prior-mean']:
            value = self.prior.mean()
        else:
            raise ValueError('`value_option` ({}) not recognized'.format(value_option))
        self.value = value

        # Set the propsal parameters
        if pl.isstr(target_acceptance_rate):
            if target_acceptance_rate in ['optimal', 'auto']:
                target_acceptance_rate = 0.44
            else:
                raise ValueError('`target_acceptance_rate` ({}) not recognized'.format(
                    target_acceptance_rate))
        elif pl.isfloat(target_acceptance_rate):
            if target_acceptance_rate < 0 or target_acceptance_rate > 1:
                raise ValueError('`target_acceptance_rate` ({}) out of range'.format(
                    target_acceptance_rate))
        else:
            raise TypeError('`target_acceptance_rate` ({}) type not recognized'.format(
                type(target_acceptance_rate)))
        self.target_acceptance_rate = target_acceptance_rate

        if pl.isstr(tune):
            if tune in ['auto']:
                tune = 50
            else:
                raise ValueError('`tune` ({}) not recognized'.format(tune))
        elif pl.isint(tune):
            if tune < 0:
                raise ValueError('`tune` ({}) must be > 0'.format(
                    tune))
        else:
            raise TypeError('`tune` ({}) type not recognized'.format(type(tune)))
        self.tune = tune

        if pl.isstr(end_tune):
            if end_tune in ['auto', 'half-burnin']:
                end_tune = int(self.G.inference.burnin/2)
            else:
                raise ValueError('`tune` ({}) not recognized'.format(end_tune))
        elif pl.isint(end_tune):
            if end_tune < 0 or end_tune > self.G.inference.burnin:
                raise ValueError('`end_tune` ({}) out of range (0, {})'.format(
                    end_tune, self.G.inference.burnin))
        else:
            raise TypeError('`end_tune` ({}) type not recognized'.format(type(end_tune)))
        self.end_tune = end_tune

        # Set the proposal variance
        if not pl.isstr(proposal_option):
            raise TypeError('`proposal_option` ({}) must be a str'.format(
                type(proposal_option)))
        elif proposal_option == 'manual':
            if not pl.isnumeric(proposal_var):
                raise TypeError('`proposal_var` ({}) must be a numeric'.format(
                    type(proposal_var)))
            if proposal_var <= 0:
                raise ValueError('`proposal_var` ({}) not proper'.format(proposal_var))
        elif proposal_option in ['auto']:
            proposal_var = (self.value ** 2)/10
        else:
            raise ValueError('`proposal_option` ({}) not recognized'.format(
                proposal_option))
        self.proposal.var.value = proposal_var
        self.proposal.low = 0
        self.proposal.high = float('inf')

    def update_var(self):
        '''Update the variance of the proposal
        '''
        if self.sample_iter == 0:
            self.temp_acceptances = 0
            self.acceptances = np.zeros(self.G.inference.n_samples, dtype=bool)
        
        elif self.sample_iter > self.end_tune:
            # Don't do any more updates
            return
        
        elif self.sample_iter % self.tune == 0:
            # Update var
            acceptance_rate = self.temp_acceptances / self.tune
            if acceptance_rate > self.target_acceptance_rate:
                self.proposal.var.value *= 1.5
            else:
                self.proposal.var.value /= 1.5
            self.temp_acceptances = 0

    def update(self):
        '''First we update the proposal (if necessary) and then we do a MH step
        '''
        self.update_var()
        proposal_std = np.sqrt(self.proposal.var.value)

        # Get the data
        xs = []
        for ridx, tidx in self.data_locs:
            xs.append(self.G[REPRNAMES.QPCR_VARIANCES].value[ridx].value[tidx])

        # Get the dof
        dof = self.G[REPRNAMES.QPCR_DOFS].value[self.l].value

        # Propose a new value for the scale
        prev_scale = self.value
        self.proposal.mean.value = self.value
        new_scale = self.proposal.sample()

        # Calculate the target distribution log likelihood
        prev_target_ll = pl.random.sics.logpdf(value=prev_scale, 
            dof=self.prior.dof.value, scale=self.prior.scale.value)
        for x in xs:
            prev_target_ll += pl.random.sics.logpdf(value=x,
                scale=prev_scale, dof=dof)
        new_target_ll = pl.random.sics.logpdf(value=new_scale, 
            dof=self.prior.dof.value, scale=self.prior.scale.value)
        for x in xs:
            new_target_ll += pl.random.sics.logpdf(value=x,
                scale=new_scale, dof=dof)

        # Normalize by the loglikelihood of the proposal
        prev_prop_ll = pl.random.truncnormal.logpdf(
            value=prev_scale, mean=new_scale, std=proposal_std,
            low=self.proposal.low, high=self.proposal.high)
        new_prop_ll = pl.random.truncnormal.logpdf(
            value=new_scale, mean=prev_scale, std=proposal_std,
            low=self.proposal.low, high=self.proposal.high)

        # Accept or reject
        r = (new_target_ll - prev_prop_ll) - \
            (prev_target_ll - new_prop_ll)
        u = np.log(pl.random.misc.fast_sample_standard_uniform())

        # print('\n\n\n{} prior_mean\n----------'.format(self.child_name))
        # print('x', x)
        # print('prev_scale', prev_scale)
        # print('prev_target_ll', prev_target_ll)
        # print('prev_prop_ll', prev_prop_ll)
        # print('new scale', new_scale)
        # print('new_target_ll', new_target_ll)
        # print('new_prop_ll', new_prop_ll)
        # print('\nr', r, u)
            
        if r >= u:
            self.acceptances[self.sample_iter] = True
            self.value = new_scale
            self.temp_acceptances += 1
        else:
            self.value = prev_scale

