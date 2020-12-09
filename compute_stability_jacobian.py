'''Compute the jacobian of the stability for each Gibb step

    ################################################
    ##
    ##
    ## CHANGE
    ##
    ################################################
    let x = -inv(A)@r
    then stability_jacobian = diag(x) @ A

Author: David Kaplan
Date: 11/30/20
MDSINE2 version: 4.0.6
'''
import mdsine2 as md2
import numpy as np
import argparse
import logging
from mdsine2.names import STRNAMES
import time

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--chain', '-c', type=str, dest='chain',
        help='This is the path of the chain for inference.')
    parser.add_argument('--outfile', '-o', type=str, dest='outfile',
        help='This is where you are saving the posterior renderings')
    parser.add_argument('--section', '-s', type=str, dest='section',
        help='Section to plot the variables of. Options: (`posterior`, ' \
            '`burnin`, `entire`)', default='posterior')
    args = parser.parse_args()
    md2.config.LoggingConfig(level=logging.INFO)

    mcmc = md2.BaseMCMC.load(args.chain)
    section = args.section

    growth = mcmc.graph[STRNAMES.GROWTH_VALUE].get_trace_from_disk(section=section)
    si = mcmc.graph[STRNAMES.SELF_INTERACTION_VALUE].get_trace_from_disk(section=section)
    interactions = mcmc.graph[STRNAMES.INTERACTIONS_OBJ].get_trace_from_disk(section=section)

    interactions[np.isnan(interactions)] = 0
    for i in range(len(mcmc.graph.data.taxas)):
        interactions[:, i, i] = -np.absolute(si[:, i])

    jacobian = np.zeros(shape=interactions.shape)

    # Get the steady state
    study = mcmc.graph.data.subjects
    steady_state = []
    for subj in study:
        M = subj.matrix()['abs']
        tidx_start = np.searchsorted(subj.times, 14)
        tidx_end = np.searchsorted(subj.times, 20)
        steady_state.append(np.mean(M[:, tidx_start:tidx_end], axis=1))
    steady_state = np.asarray(steady_state)
    steady_state = np.mean(steady_state, axis=0)
    steady_state[steady_state == 0] = 1e5
    print(steady_state.shape)
    for s in steady_state:
        print(s)

    eigan = np.zeros(shape=(interactions.shape[0], interactions.shape[1]))

    start_time = time.time()
    for gibb in range(interactions.shape[0]):
        if gibb % 5 == 0:
            logging.info('{}/{} - {}'.format(gibb, 
                interactions.shape[0], time.time()-start_time))
            start_time = time.time()
        
        r = growth[gibb].reshape(-1,1)
        A = interactions[gibb]
        # Analytical steady state
        # x_star = (-np.linalg.pinv(A) @ r).ravel()
        # diag_x_star = np.diag(x_star)
        # jacobian[gibb] = diag_x_star @ A

        # print(x_star[0])

        # # Data steady state
        # x_star = steady_state
        # diag_x_star = np.diag(x_star)

        # # Forward simulate
        # dyn = md2.model.gLVDynamicsSingleClustering(growth=r, interactions=A,
        #     sim_max=1e20) 
        # # dict ( times-> np.ndarray, 'X' -> (n_otus, n_times))
        # x = md2.integrate(dynamics=dyn, initial_conditions=steady_state.reshape(-1,1), 
        #     dt=0.01, n_days=60, subsample=True, times=np.arange(61)) 
        # x_star = x['X'][:, -1]
        # times = x['times']
        # traj = x['X']
        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        # for oidx in range(traj.shape[0]):
        #     ax.plot(times, traj[oidx, :])
        # ax.set_yscale('log')
        # plt.savefig('plt{}.pdf'.format(gibb))
        # plt.close()
        # if gibb == 100:
        #     sys.exit()
        # diag_x_star = np.diag(x_star)
        # jacobian[gibb] = np.diag((r + A@(x_star.reshape(-1,1))).ravel()) + \
        #     diag_x_star @ A

        jacobian[gibb] = np.diag(growth[gibb]) @ A
        # eigan[gibb] = np.linalg.eig(jacobian[gibb])
        
        # time.sleep(1)

    np.save(args.outfile, jacobian)
    # np.save(args.outfile, eigan)
    # import 
    # print(min_diag_xstar)


