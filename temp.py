import mdsine2 as md2
from mdsine2.names import STRNAMES
import pandas as pd
import logging
import re
import sys
import numpy as np
import matplotlib.pyplot as plt

md2.LoggingConfig(level=logging.INFO)
# study = md2.Study.load('processed_data/gibson_healthy_agg.pkl')

# for otu in study.taxas:

#     if md2.isotu(otu):
#         print(otu.name)
#         data = []
#         for asvname in otu.aggregated_taxas:
#             temp = []
#             for k in otu.aggregated_taxonomies[asvname]:
#                 temp.append(otu.aggregated_taxonomies[asvname][k])
#             data.append(temp)

#         df = pd.DataFrame(data, columns=list(otu.aggregated_taxonomies[asvname].keys()))
#         print(df)
#         print(list(otu.taxonomy.values()))

# study.taxas.generate_consensus_taxonomies()



#  rows, taxa fmt
# columns are the clusters

study = md2.Study.load('processed_data/gibson_healthy_agg.pkl')

# for taxa in study.taxas:
#     print(taxa)

for subj in study:
    for taxa in study.taxas:
        if not md2.isotu(taxa):
            continue

        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(111)
        ax = md2.visualization.aggregate_taxa_abundances(subj=subj, agg=taxa, dtype='rel', ax=ax)
        fig = plt.gcf()
        fig.tight_layout()
        # plt.savefig('output/agglomerates/{}.pdf'.format(taxa.name))
        plt.show()
