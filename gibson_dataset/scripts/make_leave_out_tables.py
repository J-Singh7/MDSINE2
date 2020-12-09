'''Make the standard leave out tables to the perturbation analysis and keystoneness
'''
import mdsine2 as md2
from mdsine2.names import STRNAMES
import argparse
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--chain', type=str, dest='chain',
        help='Location of chain')
    parser.add_argument('--output-basepath', '-o', type=str, dest='basepath',
        help='Where to save the output')
    parser.add_argument('--sep', type=str, dest='sep', default=',',
        help='separator for the leave out table')
    args = parser.parse_args()

    basepath = args.basepath
    os.makedirs(basepath, exist_ok=True)

    mcmc = md2.BaseMCMC.load(args.chain)
    study = mcmc.graph.data.subjects
    taxas = mcmc.graph.data.taxas

    # Make the taxas
    s = '\n'.join([str(i) for i in range(len(taxas))])
    fname = os.path.join(basepath, '{}-taxas.csv'.format(study.name))
    f = open(fname, 'w')
    f.write(s)
    f.close()

    # Make the clusters
    clustering = mcmc.graph[STRNAMES.CLUSTERING_OBJ]
    md2.generate_cluster_assignments_posthoc(clustering=clustering, set_as_value=True)

    ss = []
    for cluster in clustering:
        mems = [str(i) for i in cluster.members]
        ss.append(','.join(mems))
    s = '\n'.join(ss)
    fname = os.path.join(basepath, '{}-clusters.csv'.format(study.name))
    f = open(fname, 'w')
    f.write(s)
    f.close()