'''Make stacked barplots of the taxonomy in the raw data

Author: David Kaplan
Date: 11/18/20
MDSINE2 version: 4.0.6

Parameters
----------
'''
import mdsine2 as md2
import pandas as pd
import numpy as np
import sys

import matplotlib.pyplot as plt
import seaborn as sns


# Taxonomy
TAXLEVEL = 'family'
TAXLEVEL_PLURALS = {
    'genus': 'Genera', 'Genus': 'Genera',
    'family': 'Families', 'Family': 'Families',
    'order': 'Orders', 'Order': 'Orders',
    'class': 'Classes', 'Class': 'Classes',
    'phylum': 'Phyla', 'Phylum': 'Phyla',
    'kingdom': 'Kingdoms', 'Kingdom': 'Kingdoms'
}

TAXLEVEL_INTS = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
TAXLEVEL_REV_IDX = {
    'species': 0,
    'genus': 1,
    'family': 2,
    'order': 3,
    'class': 4,
    'phylum': 5,
    'kingdom': 6
}

# Aggregation of abundances below this threhold
CUTOFF_FRAC_ABUNDANCE = 0.01

# Expanded color pallete
XKCD_COLORS1 = sns.color_palette('muted', n_colors=10)
XKCD_COLORS2 = sns.color_palette("dark", n_colors=10)

XKCD_COLORS = []
for lst in [XKCD_COLORS1, XKCD_COLORS2]:
    # lst = lst[::-1]
    for c in lst:
        XKCD_COLORS.append(c)

DATA_FIGURE_COLORS = {}
XKCD_COLORS_IDX = 0

def _cnt_times(df, times, times_cnts, t2idx):
    for col in df.columns:
        if col in times:
            times_cnts[t2idx[col]] += 1
    return times_cnts

def _add_unequal_col_dataframes(df, dfother, times, times_cnts, t2idx):
    '''Add the contents of both the dataframes. This controls for the
    columns in the dataframes `df` and `dfother` being different.
    '''

    times_cnts = _cnt_times(dfother, times, times_cnts, t2idx)
    if df is None:
        return dfother, times_cnts

    cols_toadd_df = []
    cols_toadd_dfother = []
    for col in dfother.columns:
        if col not in df.columns:
            cols_toadd_df.append(col)
    for col in df.columns:
        if col not in dfother.columns:
            cols_toadd_dfother.append(col)

    df = pd.concat([df,
        pd.DataFrame(np.zeros(shape=(len(df.index), len(cols_toadd_df))), 
            index=df.index, columns=cols_toadd_df)], axis=1)
    dfother = pd.concat([dfother,
        pd.DataFrame(np.zeros(shape=(len(dfother.index), len(cols_toadd_dfother))), 
            index=dfother.index, columns=cols_toadd_dfother)], axis=1)

    return dfother.reindex(df.index) + df, times_cnts

def _get_top(df, cutoff_frac_abundance, taxlevel, taxaname_map=None):
    matrix = df.values
    abunds = np.sum(matrix, axis=1)
    namemap = {}

    a = abunds / abunds.sum()
    a = np.sort(a)[::-1]

    cutoff_num = None
    for i in range(len(a)):
        if a[i] < cutoff_frac_abundance:
            cutoff_num = i
            break
    if cutoff_num is None:
        raise ValueError('Error')
    print('Cutoff Num:', cutoff_num)

    idxs = np.argsort(abunds)[-cutoff_num:][::-1]
    dfnew = df.iloc[idxs, :]

    if taxaname_map is not None:
        indexes = df.index
        for idx in idxs:
            namemap[indexes[idx]] = taxaname_map[indexes[idx]]

    # Add everything else as 'Other'
    vals = None
    for idx in range(len(df.index)):
        if idx not in idxs:
            if vals is None:
                vals = df.values[idx, :]
            else:
                vals += df.values[idx, :]
        
    dfother = pd.DataFrame([vals], columns=df.columns, index=['{} with <{}% total abund'.format(
        TAXLEVEL_PLURALS[taxlevel], cutoff_frac_abundance*100)])
    df = dfnew.append(dfother)

    if taxaname_map is not None:
        print('\n\n')
        print('Name map')
        print('--------')
        for k,v in namemap.items():
            print()
            print(k)
            print(v)

    return df

def _make_full_df(dset):
    sys.path.append('..')
    from util import load_gibson_dataset

    subjset = load_gibson_dataset(dset=dset)
    # Make the data frame
    taxidx = TAXLEVEL_REV_IDX[TAXLEVEL]
    upper_tax = TAXLEVEL_INTS[taxidx+1]
    lower_tax = TAXLEVEL_INTS[taxidx]

    df = None
    times = []
    for subj in subjset:
        times = np.append(times, subj.times)
    times = np.sort(np.unique(times))
    t2idx = {}
    for i,t in enumerate(times):
        t2idx[t] = i
    times_cnts = np.zeros(len(times))

    for subj in subjset:
        dfnew, taxaname_map = subj.cluster_by_taxlevel(dtype='abs', lca=False, taxlevel=TAXLEVEL, 
            index_formatter='%({})s %({})s'.format(upper_tax, lower_tax), smart_unspec=False)
        df, times_cnts = _add_unequal_col_dataframes(df=df, dfother=dfnew, times=times, 
            times_cnts=times_cnts, t2idx=t2idx)

    df = df / df.sum(axis=0)

    # Only plot the OTUs that have a totol percent abundance over a threshold
    if CUTOFF_FRAC_ABUNDANCE is not None:
        df = _get_top(df, cutoff_frac_abundance=CUTOFF_FRAC_ABUNDANCE, taxlevel=TAXLEVEL)

    return df, taxaname_map

def _data_figure_rel_and_qpcr(dset, df, axqpcr, 
    axrel, axpert, axinoculum, taxaname_map, figlabelinoculum=None, figlabelqpcr=None,
    figlabelrel=None, make_legend=False, make_ylabels=True):
    '''Data summarization figure for the paper
    '''
    global DATA_FIGURE_COLORS
    global XKCD_COLORS_IDX
    subjset = md2.dataset.gibson(dset=dset)

    taxidx = TAXLEVEL_REV_IDX[TAXLEVEL]
    upper_tax = TAXLEVEL_INTS[taxidx+1]
    lower_tax = TAXLEVEL_INTS[taxidx]

    # Plot
    labels = np.asarray(list(df.index))
    labels = labels[::-1]
    matrix = df.values
    matrix = np.flipud(matrix)
    times = np.asarray(list(df.columns))

    print('Plot healthy', dset)
    print(labels)
    print(times)

    # transform the indicies for the perturbations (we're in index space not time space)
    for perturbation in subjset.perturbations:
        start_t = perturbation.start
        end_t = perturbation.end

        startidx = np.searchsorted(times, start_t)
        endidx = np.searchsorted(times, end_t)
        perturbation.start = startidx - 0.5
        perturbation.end = endidx + 0.5

    # Plot relative abundance, Create a stacked bar chart
    offset = np.zeros(matrix.shape[1])
    # colors = sns.color_palette('muted')

    for row in range(matrix.shape[0]):
        label = labels[row]
        if label in DATA_FIGURE_COLORS:
            color = DATA_FIGURE_COLORS[label]
        else:
            color = XKCD_COLORS[XKCD_COLORS_IDX]
            XKCD_COLORS_IDX += 1
            DATA_FIGURE_COLORS[label] = color

        axrel.bar(np.arange(len(times)), matrix[row,:], bottom=offset, color=color, label=label,
            width=1)
        offset = offset + matrix[row,:]

    # Set the xlabels
    locs = np.arange(0, len(times),step=10)
    ticklabels= times[locs]
    axrel.set_xticks(locs)
    axrel.set_xticklabels(ticklabels)
    axrel.yaxis.set_major_locator(plt.NullLocator())
    axrel.yaxis.set_minor_locator(plt.NullLocator())
    for tick in axrel.xaxis.get_major_ticks():
        tick.label.set_fontsize(24)

    # Plot the qpcr
    qpcr_meas = {}
    for subj in subjset:
        for t in subj.times:
            if t not in qpcr_meas:
                qpcr_meas[t] = []
            qpcr_meas[t].append(subj.qpcr[t].mean())

    # Plot the geometric mean of the subjects' qpcr
    for key in qpcr_meas:
        vals = qpcr_meas[key]
        a = 1
        for val in vals:
            a *= val
        a = a**(1/len(vals))
        qpcr_meas[key] = a

    # Sort and set
    times_qpcr = np.sort(list(qpcr_meas.keys()))
    vals = np.zeros(len(times_qpcr))
    for iii,t  in enumerate(times_qpcr):
        vals[iii] = qpcr_meas[t]
    axqpcr.plot(np.arange(len(times_qpcr)), vals, marker='o', linestyle='-', color='black')
    axqpcr.xaxis.set_major_locator(plt.NullLocator())
    axqpcr.xaxis.set_minor_locator(plt.NullLocator())
    max_qpcr_value = np.max(vals)
    # axqpcr.set_yscale('log')

    # Plot the inoculum
    # #################
    inoculum_subjset = md2.dataset.gibson(dset='inoculum')
    if dset == 'healthy':
        inoculum = inoculum_subjset['Healthy']
    else:
        inoculum = inoculum_subjset['Ulcerative Colitis']

    # print(inoculum.df()['raw'].head())
    print('Inoculum')
    df, taxaname_map_inoc = inoculum.cluster_by_taxlevel(dtype='raw', lca=False, taxlevel=TAXLEVEL,
        index_formatter='%({})s %({})s'.format(upper_tax, lower_tax), smart_unspec=False)
    df = _get_top(df, cutoff_frac_abundance=CUTOFF_FRAC_ABUNDANCE, 
        taxlevel=TAXLEVEL, taxaname_map=taxaname_map_inoc)

    matrix = df.to_numpy()
    matrix = np.flipud(matrix)
    matrix = matrix / np.sum(matrix)
    labels = np.asarray(list(df.index))
    labels = labels[::-1]

    offset = 0
    for row in range(matrix.shape[0]):
        label = labels[row]
        if label in DATA_FIGURE_COLORS:
            color = DATA_FIGURE_COLORS[label]
        else:
            color = XKCD_COLORS[XKCD_COLORS_IDX]
            XKCD_COLORS_IDX += 1
            DATA_FIGURE_COLORS[label] = color
        axinoculum.bar([0], matrix[row], bottom=[offset], label=labels, width=1, color=color)
        offset += matrix[row,0]
    axinoculum.xaxis.set_major_locator(plt.NullLocator())
    axinoculum.xaxis.set_minor_locator(plt.NullLocator())
    axinoculum.set_ylim(bottom=0, top=1)

    for tick in axinoculum.yaxis.get_major_ticks():
        tick.label.set_fontsize(24)


    # fig.subplots_adjust(left=0.07, right=0.82, hspace=0.1)
    axqpcr.set_ylabel('CFUs/g', size=25, fontweight='bold')
    # axqpcr.yaxis.set_label_coords(-0.06, 0.5)
    axinoculum.set_ylabel('Relative Abundance', size=25, fontweight='bold')
    axrel.set_xlabel('Time (d)', size=25, fontweight='bold')
    # axrel.xaxis.set_label_coords(0.5,-0.1, transform=axrel.transAxes)
    axrel.set_ylim(bottom=0,top=1)
    
    
    if dset == 'healthy':
        title = 'Healthy Cohort'
    else:
        title = 'Ulcerative Colitis Cohort'
    axqpcr.set_title(title, fontsize=30, fontweight='bold', y=1.3)
        # transform=axqpcr.transAxes)

    axpert.set_xlim(axrel.get_xlim())
    pl.visualization.shade_in_perturbations(axpert, perturbations=subjset.perturbations, textsize=22, 
        alpha=0)
    for perturbation in subjset.perturbations:
        axpert.axvline(x=perturbation.start, color='black', linestyle='--', lw=2)
        axpert.axvline(x=perturbation.end, color='black', linestyle='--', linewidth=2)

    if figlabelinoculum is not None:
        axinoculum.text(x = 0, y= 1.01,
            s=figlabelinoculum, fontsize=25, fontweight='bold',
            transform=axinoculum.transAxes)
    if figlabelqpcr is not None:
        axpert.text(x = 0, y= 1.01,
            s=figlabelqpcr, fontsize=25, fontweight='bold',
            transform=axpert.transAxes)
    if figlabelrel is not None:
        axrel.text(x = 0, y= 1.01,
            s=figlabelrel, fontsize=25, fontweight='bold',
            transform=axrel.transAxes)

    return max_qpcr_value


if __name__ == '__main__':
    global DATA_FIGURE_COLORS
    global XKCD_COLORS

    df_healthy, taxaname_map_healthy = _make_full_df('healthy')
    df_unhealthy, taxaname_map_unhealthy = _make_full_df('uc')

     # Set the colors from most to least abundant - only consider healthy
    M = df_healthy.to_numpy()
    a = np.sum(M, axis=1)
    idxs = np.argsort(a)[::-1] # reverse the indexes so it goes from largest to smallest

    taxaname_map_keys = list(taxaname_map_healthy.keys())
    print('Label, beginning')
    for idx in idxs:
        label = df_healthy.index[idx]
        color = XKCD_COLORS[XKCD_COLORS_IDX]
        XKCD_COLORS_IDX += 1
        DATA_FIGURE_COLORS[label] = color
        print(label)

    fig = plt.figure(figsize=(34,12))
    squeeze = 2
    gs = fig.add_gridspec(9,40*squeeze)

    axqpcr1 = fig.add_subplot(gs[2:4,1*squeeze:14*squeeze])
    axrel1 = fig.add_subplot(gs[4:8,1*squeeze:14*squeeze])
    axpert1 = fig.add_subplot(gs[2:8,1*squeeze:14*squeeze], facecolor='none')
    axinoculum1 = fig.add_subplot(gs[4:8,0])

    max_qpcr_value1 = _data_figure_rel_and_qpcr('healthy', df=df_healthy,
        axqpcr=axqpcr1, axrel=axrel1, axpert=axpert1,
        axinoculum=axinoculum1, make_ylabels=True,
        figlabelinoculum='D', figlabelqpcr='B', figlabelrel='E',
        make_legend=False, taxaname_map=taxaname_map_healthy)
    
