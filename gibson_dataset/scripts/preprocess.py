'''Preprocess (aggregate and filter) the Gibson dataset for Healthy cohort, 
Ulcerative Colitis cohort, inoculum, and replicate read datasets.

Author: David Kaplan
Date: 11/30/20
MDSINE2 version: 4.0.6

Methodology
-----------
1) Load the dataset
2) Aggregate the ASVs into OTUs using the aligned 16S v4 rRNA sequences in 
   `files/gibson_16S_rRNA_v4_seqs_aligned_filtered.fa` given a hamming-distance. 
   Once we agglomerate them together we set the sequences to the original sequence 
   (unaligned).
3) Calculate the consensus sequences
4) Rename the Taxas to OTUs
5) Remove selected timepoints

The file `paper_files/preprocessing/gibson_16S_rRNA_v4_seqs_aligned_filtered.fa` 
was prepared by first aligning the Taxa sequences to the reference sequeunces in the 
phylogenetic tree. Once aligned, Taxas were manually filtered out if they had poor alignment 
within the 16S rRNA v4 region. A fasta file of the Taxas removed as well as their alignments 
can be found in `paper_files/preprocessing/prefiltered_asvs.fa`. 

'''
import argparse
import pandas as pd
from Bio import SeqIO, SeqRecord, Seq
import numpy as np
import logging
import mdsine2 as md2
import os
import sys

sys.path.append('..')
import util

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--output-basepath', '-o', type=str, dest='basepath',
        help='This is where you want to save the parsed dataset.')
    parser.add_argument('--hamming-distance', '-hd', type=int, dest='hamming_distance',
        help='This is the hamming radius to aggregate Taxa sequences. If nothing ' \
            'is provided, then there will be no aggregation.', default=None)
    parser.add_argument('--rename-prefix', '-rp', type=str, dest='rename_prefix',
        help='This is the prefix you are renaming the aggregate taxas to. ' \
            'If nothing is provided, then they will not be renamed', default=None)
    parser.add_argument('--sequences', '-s', type=str, dest='sequences',
        help='This is the fasta file location of the aligned sequences for each Taxa' \
            ' that was used for placement in the phylogenetic tree. If nothing is ' \
            'provided, then do not replace them.', default=None)
    parser.add_argument('--remove-timepoints', dest='remove_timepoints', nargs='+', default=None, 
        type=float, help='Which times to remove')

    args = parser.parse_args()
    md2.config.LoggingConfig(level=logging.DEBUG)
    os.makedirs(args.basepath, exist_ok=True)

    for dset in ['healthy', 'uc', 'replicates', 'inoculum']:
        # 1) Load the dataset
        study = util.load_gibson_dataset(dset=dset, as_df=False, species_assignment='both')

        # 2) Set the sequences for each Taxa
        #    Remove all taxas that are not contained in that file
        #    Remove the gaps
        if args.sequences is not None:
            logging.info('Replacing sequences with the file {}'.format(args.sequences))
            seqs = SeqIO.to_dict(SeqIO.parse(args.sequences, format='fasta'))
            to_delete = []
            for taxa in study.taxas:
                if taxa.name not in seqs:
                    to_delete.append(taxa.name)
            for name in to_delete:
                logging.info('Deleting {} because it was not in {}'.format(
                    name, args.sequences))
            study.pop_taxas(to_delete)

            M = []
            for taxa in study.taxas:
                seq = list(str(seqs[taxa.name].seq))
                M.append(seq)
            M = np.asarray(M)
            gaps = M == '-'
            n_gaps = np.sum(gaps, axis=0)
            idxs = np.where(n_gaps == 0)[0]
            logging.info('There are {} positions where there are no gaps out of {}. Setting those ' \
                'to the sequences'.format(len(idxs), M.shape[1]))
            M = M[:, idxs]
            for i,taxa in enumerate(study.taxas):
                taxa.sequence = ''.join(M[i])

        # Aggregate with specified hamming distance
        if args.hamming_distance is not None:
            logging.info('Aggregating Taxas with a hamming distance of {}'.format(args.hamming_distance))
            study = md2.aggregate_items(subjset=study, hamming_dist=args.hamming_distance)

            # Get the maximum distance of all the OTUs
            m = -1
            for taxa in study.taxas:
                if md2.isotu(taxa):
                    for aname in taxa.aggregated_taxas:
                        for bname in taxa.aggregated_taxas:
                            if aname == bname:
                                continue
                            aseq = taxa.aggregated_seqs[aname]
                            bseq = taxa.aggregated_seqs[bname]
                            d = md2.diversity.beta.hamming(aseq, bseq)
                            if d > m:
                                m = d
            logging.info('Maximum distance within an OTU: {}'.format(m))

        # 3) compute consensus sequences
        if args.sequences is not None:
            # put original sequences in study
            orig = util.load_gibson_dataset(dset=dset, as_df=False, species_assignment='both')
            for taxa in study.taxas:
                if md2.isotu(taxa):
                    for asvname in taxa.aggregated_taxas:
                        taxa.aggregated_seqs[asvname] = orig.taxas[asvname].sequence
                else:
                    taxa.sequence = orig.taxas[taxa.name].sequence

            # Compute consensus sequences
            study.taxas.generate_consensus_seqs(threshold=0.65, noconsensus_char='N')

        # 4) Rename taxas
        if args.rename_prefix is not None:
            print('Renaming Taxas with prefix {}'.format(args.rename_prefix))
            study.taxas.rename(prefix=args.rename_prefix, zero_based_index=False)

        # 5) Remove timepoints
        if args.remove_timepoints is not None:
            if dset in ['healthy', 'uc']:
                study.pop_times(args.remove_timepoints)

        # 6) Save the study set and sequences
        study.save(os.path.join(args.basepath, 'gibson_' + dset + '_agg.pkl'))
        ret = []
        for taxa in study.taxas:
            ret.append(SeqRecord.SeqRecord(seq=Seq.Seq(taxa.sequence), id=taxa.name,
                description=''))
        SeqIO.write(ret, os.path.join(args.basepath, 'gibson_' + dset + '_agg.fa'), 'fasta-2line')
