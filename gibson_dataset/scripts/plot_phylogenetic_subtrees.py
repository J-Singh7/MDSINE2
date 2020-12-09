'''Plot the phylogenetic subtree for each taxa
'''
import ete3
from ete3 import TreeStyle
import copy
import mdsine2 as md2
from Bio import SeqIO
import argparse
import os
import treeswift
import logging
import pandas as pd
import numpy as np

if __name__ == "__main__":
    

    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--output-basepath', '-o', type=str, dest='basepath',
        help='This is where you want to save the parsed dataset.')
    parser.add_argument('--study', '-s', type=str, dest='study',
        help='Dataset that contains all of the information')
    parser.add_argument('--tree', '-t', type=str, dest='tree',
        help='Path to newick tree')
    parser.add_argument('--seq-info', type=str, dest='seq_info',
        help='Maps the sequence id of the reference sequences identifying info')
    parser.add_argument('--sep', type=str, dest='sep', default='\t')
    parser.add_argument('--family-radius-factor', type=float, dest='family_radius_factor', 
        help='How much to multiply the radius of each family', default=1.5)
    args = parser.parse_args()

    basepath = args.basepath
    os.makedirs(basepath, exist_ok=True)
    md2.LoggingConfig(level=logging.DEBUG)
    study = md2.Study.load(args.study)

    # Get the median phylogenetic distance within each family of the reference seqs
    # ------------------------------------------------------------------------------
    import treeswift
    # Make the distance matrix - this is a 2D dict
    logging.info('Making distance matrix')
    tree = treeswift.read_tree_newick(args.tree)
    M = tree.distance_matrix(leaf_labels=True)
    df_distance_matrix = pd.DataFrame(M)
    node_dict = tree.label_to_node()

    # Get the families of the reference trees
    logging.info('Read genbank file')
    df_seqs = pd.read_csv(args.seq_info, sep=args.sep, index_col=0)

    logging.info('get families of reference seqs')
    ref_families = {}
    for i, seq in enumerate(df_seqs.index):
        lineage = df_seqs['Lineage'][seq].split('; ')
        # print(lineage)
        l = len(lineage)
        if l != 7:
            continue
        family = lineage[-2].lower()
        if family not in ref_families:
            ref_families[family] = []
        ref_families[family].append(seq)

    # Get the distance of every taxa to the respective family in the reference seqs
    d = {}
    percents = []
    not_found = set([])
    for i, taxa in enumerate(study.taxas):
        if i % 100 == 0:
            logging.info('{}/{} - {}'.format(i,len(study.taxas), np.mean(percents)))
            percents = []

        if taxa.tax_is_defined('family'):
            family = taxa.taxonomy['family'].lower()
            if family not in d:
                d[family] = []

            if family not in ref_families:
                logging.debug('{} NOT IN REFERENCE TREE'.format(family))
                continue
            refs = ref_families[family]

            aaa = 0
            for ref in refs:
                try:
                    dist = tree.distance_between(node_dict[taxa.name],node_dict[ref])
                    d[family].append(dist)
                    aaa += 1
                except Exception as e:
                    not_found.add(ref)
                    # logging.debug('no worked - {}, {}'.format(taxa.name, ref))
                    continue
            percents.append(aaa/len(refs))
        else:
            print('family is not defined for', taxa.name)

    # make a text file indicating the intra family distances
    fname = os.path.join(basepath, 'family_distances.txt')
    f = open(fname, 'w')
    family_dists = {}
    for family in d:
        f.write(family + '\n')
        arr = np.asarray(d[family])
        summ = md2.summary(arr)
        family_dists[family] = summ['median']
        for k,v in summ.items():
            f.write('\t{}: {}\n'.format(k,v))

    f.write('total\n')
    arr = []
    for ele in d.values():
        arr += ele
    arr = np.asarray(arr)
    summ = md2.summary(arr)
    for k,v in summ.items():
        f.write('\t{}:{}\n'.format(k,v))

    # Set the default radius to global median
    default_radius = args.family_radius_factor * summ['median']
    f.write('default radius set to {}% of global median ({})'.format(
        100*args.family_radius_factor, default_radius))
    f.close()

    def my_layout_fn(node):
        if node.is_leaf():
            if 'OTU' in node.name:
                node.img_style["bgcolor"] = "#9db0cf"
                node.name = md2.taxaname_for_paper(taxa=study.taxas[node.name], taxas=study.taxas)
            else:
                if node.name in df_seqs.index:
                    # replace the sequence name with the species id
                    node.name = df_seqs['Species'][node.name]
                else:
                    print('Node {} not found'.format(node.name))
                

    # Make the phylogenetic subtrees for each OTU
    # -------------------------------------------
    i = 0
    names = df_distance_matrix.index
    fname = os.path.join(basepath, 'table.tsv')
    f = open(fname, 'w')
    for taxa in study.taxas:
        logging.info('\n\nLooking at {}, {}'.format(i,taxa))
        logging.info('-------------------------')
        
        f.write('{}\n'.format(taxa.name))
        tree = ete3.Tree(args.tree)
        # Get the all elements within `radius`
        

        row = df_distance_matrix[taxa.name].to_numpy()
        idxs = np.argsort(row)

        names_found = False
        mult = 1.
        title = taxa.name
        while not names_found:
            if mult == 3:
                # title += '\nTOO MANY, BREAKING'
                break
            if taxa.tax_is_defined('family'):
                family = taxa.taxonomy['family'].lower()
                if family_dists[family] < 1e-2:
                    radius = default_radius * mult
                    # title += '\nfamily defined but not ref: {} Median family distance: {:.4f}'.format(taxa.taxonomy['family'], radius)
                else:
                    radius = family_dists[family] * mult
                    # title += '\nfamily defined: Median {} distance: {:.4f}'.format(taxa.taxonomy['family'], radius)
            else:
                radius = default_radius * mult * args.family_radius_factor
                mmm = mult*100*args.family_radius_factor
                # title += '\nfamily not defined: {}% Median family distance: {:.4f}'.format(mmm, radius)

            names_to_keep = []
            for idx in idxs:
                if row[idx] > radius:
                    break
                if 'OTU' in names[idx]:
                    continue
                names_to_keep.append(names[idx])

            if len(names_to_keep) > 5:
                print(len(names_to_keep), ' found for ', taxa.name)
                names_found = True
            else:
                print('expand radius')
                # title += '\n`{}` reference seqs within radius, expanding radius by 25%'.format(len(names_to_keep))
                mult += .25

        suffix_taxa = {
            'genus': '*',
            'family': '**', 
            'order': '***', 
            'class': '****', 
            'phylum': '*****', 
            'kingdom': '******'}
        title += '\nTaxonomic Key for ASV\n'
        i = 0
        for k,v in suffix_taxa.items():
            if i == 2:
                title += '\n'
                i = 0
            title += '{}: {}'.format(v,k)
            if i == 0:
                title += ', '
            i += 1

        # print(names_to_keep)

        # Make subtree of just these names
        names_to_keep.append(taxa.name)
        tree.prune(names_to_keep, preserve_branch_length=False)

        ts = TreeStyle()
        ts.layout_fn = my_layout_fn
        ts.title.add_face(ete3.TextFace(title, fsize=15), column=1)
        fname = os.path.join(basepath, '{}.pdf'.format(taxa.name))
        tree.render(fname, tree_style=ts)
    f.close()