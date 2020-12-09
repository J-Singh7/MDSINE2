'''Figure 3 of the gibson paper
'''

import matplotlib.pyplot as plt
import seaborn as sns
import mdsine2 as md2
from mdsine2.names import STRNAMES
import numpy as np
import sys
import os
import pandas as pd

import matplotlib.image as mpimg
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
from matplotlib.colors import LogNorm

import ete3
from Bio import Phylo

def phylogenetic_heatmap_gram_split(chain_healthy, chain_uc, outfile, tree_fname, taxas):
    '''Split the phylogenetic tree into gram negative and gram positive side by side.
    Gram positive is on left, gram negative is on right.

    Parameters
    ----------
    chain_healthy, chain_uc : str
        Paths to the mcmc.pkl objects for the healthy and uc cohort for a fixed
        clustering inference, respectively
    outfile : str
        Location to save the figure
    tree_fname : str
        Location to load the newick tree that shows the placement of the OTUs
    taxas : md2.TaxaSet
        This is the MDSINE2.TaxaSet object that contains all of the OTUs. This
        is used for naming
    '''
    # Load data and get taxa names
    # ---------------------------
    os.makedirs('tmp', exist_ok=True)
    names = set([])
    for chainname in [chain_healthy, chain_uc]:
        chain = md2.BaseMCMC.load(chainname)
        for otu in chain.graph.data.taxas:
            names.add(str(otu.name))
    names = list(names)

    # Split up into gram negative and gram positive
    # ---------------------------------------------
    gramneg_taxanames = []
    grampos_taxanames = []
    for name in names:
        taxa = taxas[name]
        if md2.is_gram_negative(taxa=taxa):
            gramneg_taxanames.append(name)
        else:
            grampos_taxanames.append(name)

    # Divide up the figure into a 100x100 grid
    fig = plt.figure(figsize=(40,55))
    gs = fig.add_gridspec(100,100)

    # Make the number of rows of the heatmap proportional
    # ---------------------------------------------------
    pert_nrows = 7
    grampos_nrows = 100 - pert_nrows
    gramneg_nrows = int(grampos_nrows * len(gramneg_taxanames)/len(grampos_taxanames))
    ntwk_buffer = 0
    ntwk_row_start = 10 + gramneg_nrows + ntwk_buffer
    ntwk_nrows = int((100 - ntwk_row_start)/2)

    # Make the number of columns of the heatmap proportional to healthy and uc
    tree_ncols = 24
    tree_branch_ncols = 5
    heatmap_width = 50 - tree_ncols

    mcmc = md2.BaseMCMC.load(chain_healthy)
    healthy_nclusters = len(mcmc.graph[STRNAMES.CLUSTERING_OBJ])
    mcmc = md2.BaseMCMC.load(chain_uc)
    uc_nclusters = len(mcmc.graph[STRNAMES.CLUSTERING_OBJ])
    
    uc_ncols = int(heatmap_width*uc_nclusters/(uc_nclusters + healthy_nclusters))
    healthy_ncols = heatmap_width - uc_ncols

    print('Number of columns for UC heatmaps', uc_ncols)
    print('Number of columns for healthy heatmaps', healthy_ncols)
    print('Number of columns for phylogenetic trees', tree_ncols)
    print('Number of rows for Gram positive', grampos_nrows)
    print('Number of rows for Gram negative', gramneg_nrows)
    print('Number of rows for perturbations', pert_nrows)

    # Make the subplots
    # -----------------
    # gram positive
    ax_grampos_tree = fig.add_subplot(gs[
        pert_nrows:pert_nrows+grampos_nrows,
        0         :tree_branch_ncols])
    ax_grampos_tree_full = fig.add_subplot(gs[
        pert_nrows:pert_nrows+grampos_nrows,
        0         :tree_branch_ncols], facecolor='none')
    ax_grampos_healthy_pert = fig.add_subplot(gs[
        0         :pert_nrows,
        tree_ncols:tree_ncols+healthy_ncols])
    ax_grampos_uc_pert = fig.add_subplot(gs[
        0                       :pert_nrows,
        tree_ncols+healthy_ncols:tree_ncols+healthy_ncols+uc_ncols])
    ax_grampos_pert = fig.add_subplot(gs[
        0          :pert_nrows,
        tree_ncols:tree_ncols+healthy_ncols+uc_ncols], facecolor='none')
    ax_grampos_healthy_abund = fig.add_subplot(gs[
        pert_nrows:pert_nrows+grampos_nrows,
        tree_ncols:tree_ncols+healthy_ncols])
    ax_grampos_uc_abund = fig.add_subplot(gs[
        pert_nrows              :pert_nrows+grampos_nrows,
        tree_ncols+healthy_ncols:tree_ncols+healthy_ncols+uc_ncols])
    ax_healthy_network = fig.add_subplot(gs[
        ntwk_row_start:ntwk_row_start+ntwk_nrows,
        50            :90], zorder=-50)

    # gram negative
    ax_gramneg_tree = fig.add_subplot(gs[
        pert_nrows:pert_nrows+gramneg_nrows,
        50        :50+tree_branch_ncols])
    ax_gramneg_tree_full = fig.add_subplot(gs[
        pert_nrows :pert_nrows+grampos_nrows,
        50         :50+tree_branch_ncols], facecolor='none')
    ax_gramneg_healthy_pert = fig.add_subplot(gs[
        0            :pert_nrows,
        50+tree_ncols:50+tree_ncols+healthy_ncols])
    ax_gramneg_uc_pert = fig.add_subplot(gs[
        0                          :pert_nrows,
        50+tree_ncols+healthy_ncols:50+tree_ncols+healthy_ncols+uc_ncols])
    ax_gramneg_pert = fig.add_subplot(gs[
        0            :pert_nrows,
        50+tree_ncols:50+tree_ncols+healthy_ncols+uc_ncols], facecolor='none')
    ax_gramneg_healthy_abund = fig.add_subplot(gs[
        pert_nrows   :pert_nrows+gramneg_nrows,
        50+tree_ncols:50+tree_ncols+healthy_ncols])
    ax_gramneg_uc_abund = fig.add_subplot(gs[
        pert_nrows                  :pert_nrows+gramneg_nrows,
        50+tree_ncols+healthy_ncols:50+tree_ncols+healthy_ncols+uc_ncols])
    ax_uc_network = fig.add_subplot(gs[
        ntwk_row_start+ntwk_nrows:100,
        50                       :90], zorder=-50)

    # Plot gram positive subplots
    # ---------------------------
    ax_grampos_tree, grampos_taxaname_order = _make_phylogenetic_tree(
        tree_fname=tree_fname, names=grampos_taxanames, taxas=taxas, fig=fig,
        ax=ax_grampos_tree, figlabel='E', figlabelax=ax_grampos_tree_full)
    
    ax_grampos_healthy_abund, grampos_healthy_colorder = _make_cluster_membership_heatmap(
        chainname=chain_healthy, ax=ax_grampos_healthy_abund, order=grampos_taxaname_order,
        binary=False, fig=fig, make_colorbar=False)

    ax_grampos_uc_abund, grampos_uc_colorder = _make_cluster_membership_heatmap(
        chainname=chain_uc, ax=ax_grampos_uc_abund, order=grampos_taxaname_order,
        binary=False, fig=fig, make_colorbar=False)

    ax_grampos_healthy_pert = _make_perturbation_heatmap(chainname=chain_healthy,
        min_bayes_factor=10, ax=ax_grampos_healthy_pert, colorder=grampos_healthy_colorder,
        fig=fig, make_colorbar=False, figlabel='A', render_labels=True)
    ax_grampos_uc_pert = _make_perturbation_heatmap(chainname=chain_uc,
        min_bayes_factor=10, ax=ax_grampos_uc_pert, colorder=grampos_uc_colorder,
        fig=fig, make_colorbar=False, figlabel='B', render_labels=False)

    ax_grampos_tree = _remove_border(ax_grampos_tree)
    ax_grampos_pert = _remove_border(ax_grampos_pert)
    ax_grampos_pert.set_title('Gram + Bacteria', fontsize=50, pad=80)

    # Plot gram negative subplots
    # ---------------------------
    ax_gramneg_tree, gramneg_taxaname_order = _make_phylogenetic_tree(
        tree_fname=tree_fname, names=gramneg_taxanames, taxas=taxas, fig=fig,
        ax=ax_gramneg_tree, figlabel='F', figlabelax=ax_gramneg_tree_full)
    
    ax_gramneg_healthy_abund, gramneg_healthy_colorder = _make_cluster_membership_heatmap(
        chainname=chain_healthy, ax=ax_gramneg_healthy_abund, order=gramneg_taxaname_order,
        binary=False, fig=fig, make_colorbar=False)

    ax_gramneg_uc_abund, gramneg_uc_colorder = _make_cluster_membership_heatmap(
        chainname=chain_uc, ax=ax_gramneg_uc_abund, order=gramneg_taxaname_order,
        binary=False, fig=fig, make_colorbar=True)

    ax_gramneg_healthy_pert = _make_perturbation_heatmap(chainname=chain_healthy,
        min_bayes_factor=10, ax=ax_gramneg_healthy_pert, colorder=gramneg_healthy_colorder,
        fig=fig, make_colorbar=False, figlabel='C', render_labels=True)
    ax_gramneg_uc_pert = _make_perturbation_heatmap(chainname=chain_uc,
        min_bayes_factor=10, ax=ax_gramneg_uc_pert, colorder=gramneg_uc_colorder,
        fig=fig, make_colorbar=True, figlabel='D', render_labels=False)

    ax_gramneg_tree = _remove_border(ax_gramneg_tree)
    ax_gramneg_pert = _remove_border(ax_gramneg_pert)
    ax_gramneg_tree_full = _remove_border(ax_gramneg_tree_full)
    ax_grampos_tree_full = _remove_border(ax_grampos_tree_full)
    
    ax_gramneg_pert.set_title('Gram - Bacteria', fontsize=50, pad=80)

    # # Render the networks below the gram negative
    # # -------------------------------------------
    # # healthy
    # ax_healthy_network = _remove_border(ax_healthy_network)
    # arr_man = mpimg.imread('figures/hairball_network_healthy.jpg')
    # imagebox = OffsetImage(arr_man, zoom=0.325)
    # ab = AnnotationBbox(imagebox, (0.5, 0.475), pad=0, box_alignment=(0.5,0.5))
    # ax_healthy_network.add_artist(ab)
    # ax_healthy_network.text(x=0.3, y=0.9, s='G', fontsize=50, fontweight='bold',
    #     transform=ax_healthy_network.transAxes, zorder=50)

    # # uc
    # ax_uc_network = _remove_border(ax_uc_network)
    # arr_man = mpimg.imread('figures/hairball_network_uc.jpg')
    # imagebox = OffsetImage(arr_man, zoom=0.325)
    # ab = AnnotationBbox(imagebox, (0.5, 0.475), pad=0, box_alignment=(0.5,0.5))
    # ax_uc_network.add_artist(ab)
    # ax_uc_network.text(x=0.3, y=0.9, s='H', fontsize=50, fontweight='bold',
    #     transform=ax_uc_network.transAxes, zorder=50)
    
    # Add in the taxonomic key
    suffix_taxa = {'genus': '*',
        'family': '**', 'order': '***', 'class': '****', 'phylum': '*****', 'kingdom': '******'}
    text = '$\\bf{Taxonomy} \\bf{Key}$\n'
    for taxa in suffix_taxa:
        text += '{} - {}\n'.format(suffix_taxa[taxa], taxa)
    fig.text(x=0.875, y=0.45, s=text, fontsize=35)


    fig.subplots_adjust(wspace=0.2, left=0.015, right=0.985, hspace=0.05,
        top=.96, bottom=0.02)
    plt.savefig(outfile)
    plt.close()

def _remove_border(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_minor_locator(plt.NullLocator())
    ax.set_xlabel('')
    ax.set_ylabel('')
    return ax

def _make_phylogenetic_tree(tree_fname, names, taxas, ax, fig, figlabel=None, figlabelax=None):

    tree = ete3.Tree(tree_fname)
    tree.prune(names, True)
    tree.write(outfile='tmp/temp.nhx')
    fontsize=22

    taxonomies = ['family', 'order', 'class', 'phylum', 'kingdom']
    suffix_taxa = {'genus': '*',
        'family': '**', 'order': '***', 'class': '****', 'phylum': '*****', 'kingdom': '******'}
    extra_taxa_added = set([])

    tree = Phylo.read('tmp/temp.nhx', format='newick')
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False)
    taxa_order = []
    for text in ax.texts:
        taxa_order.append(text._text)
        # Substitute the name of the taxa with the species/genus if possible
        taxaname = str(text._text).replace(' ','')
        taxa = taxas[taxaname]
        suffix = '' # for defining taxonomic level outside genus
        if taxa.tax_is_defined('genus'):
            taxaname = ' ' + taxa.taxonomy['genus']
            if taxa.tax_is_defined('species'):
                spec = taxa.taxonomy['species']
                l = spec.split('/')
                if len(l) < 3:
                    spec = '/'.join(l)
                    taxaname = taxaname + ' {}'.format(spec)
                elif len(l) >= 3:
                    spec = '/'.join(l[:2])
                    taxaname = taxaname + ' {}'.format(spec)
            else:
                suffix = suffix_taxa['genus']
        else:
            found = False
            for taxalevel in taxonomies:
                if found:
                    break
                if taxa.tax_is_defined(taxalevel):
                    found = True
                    taxaname = ' ' + taxa.taxonomy[taxalevel]
                    suffix = suffix_taxa[taxalevel]
                    extra_taxa_added.add(taxalevel)

            if not found:
                taxaname = '#'*40

        taxaname += ' ' + taxa.name
        taxaname = ' ' + suffix + taxaname
        text._text = taxaname
        text._text = text._text + '- ' * 55
        text.set_fontsize(fontsize)

    if figlabel is not None:
        figlabelax.text(x=0.15, y=1.0, s=figlabel, fontsize=50, fontweight='bold',
            transform=ax.transAxes)

    # # Make the taxnonmic key on the right hand side
    # text = '$\\bf{Taxonomy} \\bf{Key}$\n'
    # for taxa in suffix_taxa:
    #     text += '{} - {}\n'.format(suffix_taxa[taxa], taxa)
    # if side_by_side is not None:
    #     if not side_by_side:
    #         fig.text(0.1, 0.875, text, fontsize=18)
    #     else:
    #         fig.text(0.9, 0.3, text, fontsize=18)
    return ax, taxa_order

def _make_perturbation_heatmap(chainname, min_bayes_factor, ax, colorder, fig, make_colorbar=True, figlabel=None,
    render_labels=True):
    chain = md2.BaseMCMC.load(chainname)
    subjset = chain.graph.data.subjects
    clustering = chain.graph[STRNAMES.CLUSTERING_OBJ]

    matrix = np.zeros(shape=(len(chain.graph.perturbations), len(clustering)))

    index = []
    for pidx, perturbation in enumerate(chain.graph.perturbations):
        index.append(perturbation.name)

        bayes_factor = md2.generate_perturbation_bayes_factors_posthoc(mcmc=chain, 
            perturbation=perturbation, section='posterior')
        values = md2.summary(perturbation, section='posterior', only=['mean'])['mean']

        for cidx in range(len(clustering)):
            aidx = list(clustering.clusters[clustering.order[cidx]].members)[0]
            if bayes_factor[aidx] >= min_bayes_factor:
                matrix[pidx, cidx] = values[aidx]

    df = pd.DataFrame(matrix, columns=_make_names_of_clusters(len(clustering)), 
        index=index)

    # with open('perturbation_heatmap.pkl', 'wb') as handle:
    #     pickle.dump(df, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # print(df.head())

    df = df[colorder]

    max_val = np.max(np.absolute(df.values))
    max_val = 5
    im = ax.imshow(df.values, cmap='bwr_r', aspect='auto', vmin=-max_val, 
        vmax=max_val)
    ax.yaxis.set_major_locator(plt.NullLocator())
    ax.xaxis.set_major_locator(plt.NullLocator())
    
    # Make grid
    ax.set_xticks(np.arange(0.5, len(df.columns), 1), minor=True)
    ax.set_yticks(np.arange(0.5, len(subjset.perturbations), 1), minor=True)
    ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.1)

    ax.tick_params(axis='both', which='minor', left=False, bottom=False)
    
    if make_colorbar:
        cbaxes = fig.add_axes([0.92, 0.3, 0.02, 0.1]) # left, bottom, width, height
        cbar = plt.colorbar(im, cax=cbaxes, orientation='vertical', ticks=[-5,-2.5,0,2.5,5])
        cbar.ax.set_yticklabels(['<-5', '-2.5', '0', '2.5', '>5'], fontsize=30)
        cbar.ax.set_title('Perturbation\nEffect', fontsize=35, fontweight='bold')

    if render_labels:
        ax.set_yticks(np.arange(len(subjset.perturbations)), minor=False)
        ax.set_yticklabels(list(df.index))
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_rotation(0)
            tick.label.set_fontsize(30)
    else:
        ax.set_yticks([])

    if figlabel is not None:
        ax.text(x=0.0, y=1.05, s=figlabel, fontsize=50, fontweight='bold',
            transform=ax.transAxes)

    return ax

def _make_cluster_membership_heatmap(chainname, ax, order, binary, fig, make_colorbar=True, figlabel=None):
    '''Make the heatmap of the cluster membership
    If `binary` is True, then we just have a binary vector of the membership of
    the taxa in the cluster. If `binary` is False, then the coloring is the average relative
    abundance of the ASV and the colorbar is on a log scale.
    '''
    chain = md2.BaseMCMC.load(chainname)
    subjset = chain.graph.data.subjects
    clusters = chain.graph[STRNAMES.CLUSTERING_OBJ].tolistoflists()

    # If binary is False, then get the relative abundances of the ASVS
    rel_abund = None
    if not binary:
        rel_abund = np.zeros(len(subjset.taxas))
        for subj in subjset:
            M = subj.matrix()['rel']
            start_idx = np.searchsorted(subj.times, 7)
            end_idx = np.searchsorted(subj.times, 20)

            rel_abund += np.mean(M[:,start_idx:end_idx], axis=1)
        rel_abund /= len(subjset)
        max_rel = np.max(rel_abund)
        min_rel = np.min(rel_abund)

    # print('here')
    # print(min_rel)
    if min_rel < 0.001:
        min_rel = 0.001
    rel_abund[rel_abund == 0] = min_rel
    # print(max_rel)

    matrix = np.zeros(shape=(len(subjset.taxas), len(clusters)))
    for cidx, cluster in enumerate(clusters):
        for oidx in cluster:
            if binary:
                matrix[oidx, cidx] = 1
            else:
                matrix[oidx, cidx] = rel_abund[oidx]

    index = [taxa.name for taxa in subjset.taxas]

    iii = 0
    for taxanew in order:
        if taxanew not in index:
            index.append(taxanew)
            iii += 1

    # Add in nan rows in places that order is not in index
    matrix = np.vstack((matrix, np.zeros(shape=(iii, matrix.shape[1]))*np.nan))
    order = [str(a).replace(' ', '') for a in order]

    df = pd.DataFrame(matrix,
        columns=_make_names_of_clusters(len(clusters)),
        index=index)

    
    # indices = {}
    # for name in df.index:
    #     indices[name] = name.replace('OTU', 'ASV')
    # df1 = df.rename(index=indices)
    
    # with open('cluster_heatmap.pkl', 'wb') as handle:
    #     pickle.dump(df1, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # print(df1.head())
    
    # Reindex the dataframe based on the order of the phylogenetic tree
    df = df.reindex(order)

    # # Reindex the columns based off of a dendrogram
    # values = df.values
    # Z = linkage(values.T, 'ward', optimal_ordering=True)
    # n_cols = df.shape[1]
    # cols = []
    # for row in range(Z.shape[0]):
    #     z1 = Z[row,0]
    #     z2 = Z[row,1]
    #     if z1 < n_cols:
    #         cols.append(z1)
    #     if z2 < n_cols:
    #         cols.append(z2)
    # cols = np.asarray(cols)[::-1]
    cols = np.arange(len(df.columns))

    # print(Z)
    # print(cols)
    colnames = df.columns
    newcolnames = []
    for idx in cols:
        newcolnames.append(colnames[idx])
    df = df[newcolnames]

    if not binary:
        kwargs = {'norm': LogNorm(vmin=min_rel, vmax=max_rel)}
    else:
        kwargs = {}

    cmap = sns.cubehelix_palette(n_colors=100, as_cmap=True, start=2, rot=0, dark=0, light=0.5)
    cmap.set_bad(color='silver')
    cmap.set_under(color='white')
    vals = df.values
    vals += min_rel - 1e-10
    im = ax.imshow(vals, cmap=cmap, aspect='auto', **kwargs)
    ax.yaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())

    xticklabels = np.arange(1, len(df.columns)+1, step=2, dtype=int)
    ax.set_xticks(ticks=np.arange(len(df.columns), step=2))
    ax.set_xticklabels(labels=xticklabels)
    for i, tick in enumerate(ax.xaxis.get_major_ticks()):
        tick.label.set_rotation(90)
        tick.label.set_fontsize(25)
        # tick.label.set_horizontalalignment('right')
    
    # Make grid
    ax.set_xticks(np.arange(0.5, len(df.columns), 1), minor=True)
    ax.set_yticks(np.arange(0.5, len(df.index), 1), minor=True)
    ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.1)

    ax.tick_params(axis='both', which='minor', left=False, bottom=False)

    if make_colorbar:
        cbaxes = fig.add_axes([0.92, 0.1, 0.02, 0.1]) # left, bottom, width, height
        cbar = plt.colorbar(im, cax=cbaxes, orientation='vertical')
        cbar.ax.set_title('Relative\nAbundance', fontsize=30, fontweight='bold')
        cbar.ax.tick_params(labelsize=40)

    return ax, newcolnames

def _make_names_of_clusters(n_clusters):
    '''Standardize the way to name the clusters so that they are all the same between
    eachother
    '''
    return ['Cluster {}'.format(i) for i in range(n_clusters)]

if __name__ == '__main__':

    chain_healthy = '../output/mdsine2/fixed_clustering/healthy/mcmc.pkl'
    chain_uc = '../output/mdsine2/fixed_clustering/uc/mcmc.pkl'
    outfile = 'tmp/figure_3.pdf'
    tree_fname = '../gibson_dataset/files/phylogenetic_placement_OTUs/phylogenetic_tree_only_query.nhx'
    study = md2.Study.load('../processed_data/gibson_healthy_agg_taxa.pkl')
    taxas = study.taxas

    phylogenetic_heatmap_gram_split(chain_healthy=chain_healthy, chain_uc=chain_uc, 
        outfile=outfile,
        tree_fname=tree_fname, taxas=taxas)