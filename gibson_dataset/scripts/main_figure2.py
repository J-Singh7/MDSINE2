#make figure 2
'''
python main_figure2.py -file1 "../../processed_data/gibson_healthy_agg.pkl" \
       -file2 "../../processed_data/gibson_uc_agg.pkl" \
       -file3 "../../processed_data/gibson_inoculum_agg.pkl"
'''

import mdsine2 as md2
import pandas as pd
import numpy as np
import argparse

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as plticker
from matplotlib.ticker import ScalarFormatter, LogFormatter, LogFormatterSciNotation, FixedLocator
import matplotlib.patches as patches
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
# from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import TransformedBbox, BboxPatch, BboxConnector
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.legend as mlegend
from matplotlib.colors import LogNorm
import matplotlib.image as mpimg
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
from matplotlib.gridspec import GridSpec

TAXLEVEL = "family"
TAXLEVEL_PLURALS = {'genus': 'Genera', 'Genus': 'Genera', 'family': 'Families',
                'Family': 'Families', 'order': 'Orders', 'Order': 'Orders',
                'class': 'Classes', 'Class': 'Classes', 'phylum': 'Phyla',
                'Phylum': 'Phyla', 'kingdom': 'Kingdoms', 'Kingdom': 'Kingdoms'}

PERTURBATION_COLOR = "orange"
TAXLEVEL_INTS = ["species", "genus", "family", "order", "class", "phylum",
                    "kingdom"]
TAXLEVEL_REV_IDX = {"species" : 0, "genus" : 1, "family" : 2, "order" : 3,                      "class" : 4, "phylum" : 5, "kingdom" : 6}

#Aggregation of abundances below this threhold
CUTOFF_FRAC_ABUNDANCE = 0.01

def parse_args():

    parser = argparse.ArgumentParser(description = "files needed for making"\
    "main figure 2")
    parser.add_argument("-file1", "--healthy_pkl", required = "True",
        help = "pickled pl.base.Study file for healthy subjects")
    parser.add_argument("-file2", "--uc_pkl", required = "True",
        help = "pickled pl.base.Study file for UC subjects")
    parser.add_argument("-file3", "--inoc_pkl", required = "True",
        help = "pickled pl.base.Study file for inoculum")

    return parser.parse_args()


def _cnt_times(df, times, times_cnts, t2idx):
    """counts the number of times data at a given point were collected"""

    for col in df.columns:
        if col in times:
            times_cnts[t2idx[col]] += 1
    #print(times_cnts)

    return times_cnts

def _add_unequal_col_dataframes(df, dfother, times, times_cnts, t2idx):
    '''
    Add the contents of both the dataframes. This controls for the
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
    """
       selects the data associated with taxa (at taxlevel) whose abundace is
       greater than the cutoff_frac_abundance
    """
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
    #Number of taxa whose abundance is greater than the threshold
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

    #if taxaname_map is not None:
    #    print('\n\n')
    #    print('Name map')
    #    print('--------')
    #    for k,v in namemap.items():
    #        print()
    #        print(k)
    #        print(v)

    return df

def get_df(subjset):
    """
       @parameters
       subjset : (pl.Subject)
    """
    taxidx = TAXLEVEL_REV_IDX[TAXLEVEL]
    upper_tax = TAXLEVEL_INTS[taxidx+1]
    lower_tax = TAXLEVEL_INTS[taxidx]

    df = None
    times = []
    for subj in subjset:
        times = np.append(times, subj.times)

    times = np.sort(np.unique(times))#the times at which samples were collected
    t2idx = {}
    for i,t in enumerate(times):
        t2idx[t] = i
    times_cnts = np.zeros(len(times)) #the times at which samples were taken

    #update the data frame for each subject
    for subj in subjset:
        dfnew, taxaname_map = subj.cluster_by_taxlevel(dtype='abs',
        taxlevel=TAXLEVEL, index_formatter='%({})s %({})s'.format(upper_tax,
        lower_tax), smart_unspec=False)

        df, times_cnts = _add_unequal_col_dataframes(df=df, dfother=dfnew,
             times=times, times_cnts=times_cnts, t2idx=t2idx)

    df = df / df.sum(axis=0)

    # Only plot the OTUs that have a totol percent abundance over a threshold
    if CUTOFF_FRAC_ABUNDANCE is not None:
        df = _get_top(df, cutoff_frac_abundance=CUTOFF_FRAC_ABUNDANCE,
              taxlevel=TAXLEVEL)

    return df, taxaname_map

def checker(df1, df2):

    sum1 = np.sum(df1.to_numpy(), axis = 0)
    sum2 = np.sum(df2.to_numpy(), axis = 0)

def set_colors(df, color_idx, color_taxa_dict, color_set):
    """choose the color to represent each taxa level"""

    M = df.to_numpy()
    a = np.sum(M, axis = 1)
    idxs = np.argsort(a)[::-1]
    for idx in idxs:
        #print("color:", color_idx)
        label = df.index[idx]
        color = color_set[color_idx]
        color_idx += 1
        color_taxa_dict[label] = color

    return color_idx

def add_perturbation_label(ax, perturbations, subj, times, textcolor='black',
     textsize=None, alpha=0.25, label=True):

    if md2.issubject(subj):
        subj = subj.name
    if not md2.isstr(subj):
        raise ValueError('`Cannot recognize {}'.format(subj))
    if perturbations is None or len(perturbations) == 0:
        return ax

    pert_locs = []
    pert_names = []
    for pidx, perturbation in enumerate(perturbations):

        if subj not in perturbation.starts or subj not in perturbation.ends:
            continue
        ax.axvspan(
            xmin = perturbation.starts[subj],
            xmax = perturbation.ends[subj],
            facecolor = PERTURBATION_COLOR,
            alpha=alpha, zorder=-10000)
        pert_locs.append((times.index(perturbation.ends[subj]) +
        times.index(perturbation.starts[subj])) / 2)
        name = perturbation.name
        if name is None:
            name = 'pert{}'.format(pidx)
        pert_names.append(name)

    if label:
        # Set the names on the top x-axis
        ax2 = ax.twiny()

        # # Set the visibility of the twin axis to see through
        ax2.spines['top'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.xaxis.set_major_locator(plt.NullLocator())
        ax2.xaxis.set_minor_locator(plt.NullLocator())
        ax2.yaxis.set_major_locator(plt.NullLocator())
        ax2.yaxis.set_minor_locator(plt.NullLocator())

        left,right = ax.get_xlim()
        ax2.set_xlim(ax.get_xlim())
        pl = []
        pn = []
        for idx, loc in enumerate(pert_locs):
            if loc > left and loc < right:
                pl.append(loc)
                pn.append(pert_names[idx])
        ax2.set_xticks(pl)
        ax2.set_xticklabels(pn)
        ax2.tick_params('x', which='both', length=0, colors=textcolor,
            labelsize=textsize)

    return ax

def plot_rel_and_qpcr(subjset, subjset_inoc, df, dset_type, axqpcr, axrel, axpert,
    axinoculum, taxaname_map, color_taxa_dict, color_index, color_set,
    figlabelinoculum = None, figlabelqpcr = None, figlabelrel = None,
    make_legend = False, make_ylabels = True, labels_order = None,
    inoc_order = None):
    """
    plots the relative abundance and perturbation
    """

    taxidx = TAXLEVEL_REV_IDX[TAXLEVEL]
    upper_tax = TAXLEVEL_INTS[taxidx+1]
    lower_tax = TAXLEVEL_INTS[taxidx]

    final_labels = []
    labels = None
    if labels_order is None:
        labels = np.asarray(list(df.index))
        labels = labels[::-1]
    else:
        labels = np.asarray(list(df.index))
        labels = labels[::-1]
        for lab in labels_order:
            if lab in labels:
                final_labels.append(lab)

        for lab in labels:
            if lab not in final_labels:
                final_labels.append(lab)
        labels = final_labels
    matrix = df.values
    matrix = np.flipud(matrix)
    times = np.asarray(list(df.columns))
    #print("shapes", matrix.shape)
    if labels_order is not None:
        new_df = df.reindex(final_labels[::-1])
        matrix = new_df.values
        matrix = np.flipud(matrix)
        times = np.asarray(list(new_df.columns))
        #print("shapes", matrix.shape)

    # Plot relative abundance, Create a stacked bar chart
    offset = np.zeros(matrix.shape[1])
    for row in range(matrix.shape[0]):
        label = labels[row]
        if label in color_taxa_dict:
            color = color_taxa_dict[label]
        else:
            color = color_set[color_index]
            color_index += 1
            color_taxa_dict[label] = color

        axrel.bar(np.arange(len(times)), matrix[row,:], bottom=offset,
        color=color, label=label, width=1, linewidth = 1)
        offset = offset + matrix[row,:]

    #set the xlabels
    locs = np.arange(0, len(times), step = 10)
    ticklabels = times[locs]
    axrel.set_xticks(locs)
    axrel.set_xticklabels(ticklabels)
    axrel.yaxis.set_major_locator(plt.NullLocator())
    axrel.yaxis.set_minor_locator(plt.NullLocator())
    for tick in axrel.xaxis.get_major_ticks():
        tick.label.set_fontsize(24)

    #plot the absolute abundance
    qpcr_meas = {}
    for subj in subjset:
        for t in subj.times:
            if t not in qpcr_meas:
                qpcr_meas[t] = []
            qpcr_meas[t].append(subj.qpcr[t].mean())
    #geometric mean of qpcr
    for key in qpcr_meas:
        vals = qpcr_meas[key]
        a = 1
        for val in vals:
            a *= val
        a = a ** (1 / len(vals))
        qpcr_meas[key] = a
    times_qpcr = np.sort(list(qpcr_meas.keys()))
    vals = np.zeros(len(times_qpcr))
    for iii, t in enumerate(times_qpcr):
        vals[iii] = qpcr_meas[t]

    axqpcr.plot(np.arange(len(times_qpcr)), vals, marker='o', linestyle='-',
           color='black')
    axqpcr.xaxis.set_major_locator(plt.NullLocator())
    axqpcr.xaxis.set_minor_locator(plt.NullLocator())
    max_qpcr_value = np.max(vals)

    inoc = None
    if dset_type == "healthy":
        inoc = subjset_inoc["Healthy"]
    else:
        inoc = subjset_inoc["Ulcerative Colitis"]

    print("Adding inoculum")
    df_inoc, taxa_map_inoc = inoc.cluster_by_taxlevel(dtype='raw',
            taxlevel=TAXLEVEL, index_formatter='%({})s %({})s'.format(upper_tax,
            lower_tax), smart_unspec=False)
    df_inoc = _get_top(df_inoc, cutoff_frac_abundance=CUTOFF_FRAC_ABUNDANCE,
        taxlevel=TAXLEVEL, taxaname_map = taxa_map_inoc)

    matrix_inoc = df_inoc.to_numpy()
    matrix_inoc = np.flipud(matrix_inoc)
    matrix_inoc = matrix_inoc / np.sum(matrix_inoc)
    labels_inoc = np.asarray(list(df_inoc.index))
    labels_inoc = labels_inoc[::-1]
    final_inoc_labels = []
    if inoc_order is not None:
        for lab in inoc_order:
            if lab in labels_inoc:
                final_inoc_labels.append(lab)
        for lab in labels_inoc:
            if lab not in final_inoc_labels:
                final_inoc_labels.append(lab)
        labels_inoc = final_inoc_labels
    if inoc_order is not None:
        new_df_inoc = df_inoc.reindex(final_inoc_labels[::-1])
        matrix_inoc = new_df_inoc.values
        matrix_inoc = np.flipud(matrix_inoc)
        matrix_inoc = matrix_inoc / np.sum(matrix_inoc)
    #plot the inoclum
    offset_inoc = 0
    for row in range(matrix_inoc.shape[0]):
        label = labels_inoc[row]
        if label in color_taxa_dict:
            color = color_taxa_dict[label]
        else:
            color = color_set[color_index]
            color_index += 1
            color_taxa_dict[label] = color
        axinoculum.bar([0], matrix_inoc[row], bottom=[offset_inoc],
        label = labels_inoc, width=1, color=color)
        offset_inoc += matrix_inoc[row,0]

    axinoculum.xaxis.set_major_locator(plt.NullLocator())
    axinoculum.xaxis.set_minor_locator(plt.NullLocator())
    axinoculum.set_ylim(bottom=0, top=1)

    for tick in axinoculum.yaxis.get_major_ticks():
        tick.label.set_fontsize(24)
    axqpcr.set_ylabel('CFUs/g', size=  24, fontweight='bold')
    # axqpcr.yaxis.set_label_coords(-0.06, 0.5)
    axinoculum.set_ylabel('Relative Abundance', size = 24, fontweight='bold')
    axrel.set_xlabel('Time (d)', size=  24, fontweight='bold')
    # axrel.xaxis.set_label_coords(0.5,-0.1, transform=axrel.transAxes)
    axrel.set_ylim(bottom=0, top=1)

    if dset_type == 'healthy':
        title = 'Healthy Cohort'
    else:
        title = 'Ulcerative Colitis Cohort'
    axqpcr.set_title(title, fontsize = 30, fontweight='bold', y=1.3)

    #plot perturbation
    axpert.set_xlim(axrel.get_xlim())
    subj_ = ""
    for sub in subjset:
        subj_ = sub
        break
    times_li = list(times)
    axpert.set_xlim(axrel.get_xlim())
    axpert = add_perturbation_label(axpert, subjset.perturbations, subj_, list(times),
        textsize = 22, alpha=0)

    for perturbation in subjset.perturbations:
        axpert.axvline(x = times_li.index(perturbation.starts[subj_.name]),
                       color='black', linestyle='--', lw=2)
        axpert.axvline(x = times_li.index(perturbation.ends[subj_.name]),
                      color='black', linestyle='--', linewidth=2)

    if figlabelinoculum is not None:
        print("label", figlabelinoculum)
        axinoculum.text(0, y = 1.01, s = figlabelinoculum, fontsize = 25,
                  fontweight = "bold", transform = axinoculum.transAxes)
    if figlabelqpcr is not None:
        print("label", figlabelqpcr)
        axpert.text(0, y = 1.01, s = figlabelqpcr, fontsize = 25,
                  fontweight = "bold", transform = axpert.transAxes)
    if figlabelrel is not None:
        print("label", figlabelrel)
        axrel.text(0, y = 1.01, s = figlabelrel, fontsize = 25,
                  fontweight = "bold", transform = axrel.transAxes)

    return max_qpcr_value, color_index, labels, labels_inoc

def add_expt_figure(ax, subjset, figlabel):

    times = []
    for subj in subjset:
        times = np.append(times, subj.times)
    times = np.sort(np.unique(times))

    y = [-0.1 for t in times]
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_minor_locator(plt.NullLocator())

    markerline, stemlines, baseline = ax.stem(times, y, linefmt='none')
    baseline.set_color('black')
    markerline.set_color('black')
    markerline.set_markersize(5)

    x = np.arange(0, np.max(times), step=10)
    labels = ['Day {}'.format(int(t)) for t in x]

    for ylim in [0.15, -0.35]:
        y = [ylim for t in x]
        markerline, stemlines, baseline = ax.stem(x, y)
        #stemlines = [stemline.set_color('black') for stemline in stemlines]
        stemlines.set_color("black")
        baseline.set_color('black')
        markerline.set_color('none')
    for i in range(len(labels)):
        label = labels[i]
        xpos = x[i]
        ax.text(xpos, 0.-.6, label, horizontalalignment='center', fontsize=24)
    x = np.arange(0,np.max(times),2)
    for ylim in [0.07, -0.07]:
        y = [ylim for t in x]
        markerline, stemlines, baseline = ax.stem(x, y)
        #stemlines = [stemline.set_color('black') for stemline in stemlines]
        stemlines.set_color("black")
        baseline.set_color('black')
        markerline.set_color('none')
    subj_ = None
    for subj in subjset:
        subj_ = subj
        break

    times_li = list(times)

    for perturbation in subjset.perturbations:
        name = perturbation.name
        x = (perturbation.ends[subj_.name] + perturbation.starts[subj_.name])/2
        print(x)
        ax.text(x, 0.25, name.capitalize(), horizontalalignment='center',
            fontsize = 24)
        starts = np.asarray([perturbation.starts[subj_.name]])
        ends = np.asarray([perturbation.ends[subj_.name]])
        ax.barh(y=[0 for i in range(len(starts))], width=ends-starts, height=0.1,
        left=starts, color='darkgrey')
    ax.text(0, 0.25, 'Colonization', horizontalalignment = 'center', fontsize=24)
    if figlabel is not None:
        ax.text(x=-0.08, y = 3, s=figlabel, fontsize=25, fontweight='bold',
            transform = ax.transAxes)
    xpos = np.max(times)* 1.05
    y = 0.05
    ax.scatter([xpos], [y], c = 'black', s=25)
    ax.text(xpos+1, y, 'Fecal Sample Collection', horizontalalignment='left',
    fontsize = 24, verticalalignment='center')


def main():

    XKCD_COLORS1 = sns.color_palette('muted', n_colors=10)
    XKCD_COLORS2 = sns.color_palette("dark", n_colors=10)

    #get more colors
    XKCD_COLORS = []
    for lst in [XKCD_COLORS1, XKCD_COLORS2]:
        # lst = lst[::-1]
        for c in lst:
            XKCD_COLORS.append(c)

    DATA_FIGURE_COLORS = {}
    XKCD_COLORS_IDX = 0
    SAVEPATH = "output_figures/"

    args = parse_args()
    subjset_healthy = md2.Study.load(args.healthy_pkl)
    subjset_uc = md2.Study.load(args.uc_pkl)
    subjset_inoc = md2.Study.load(args.inoc_pkl)

    df_healthy, taxa_map_healthy = get_df(subjset_healthy)
    df_uc, taxa_map_uc = get_df(subjset_uc)
    XKCD_COLORS_IDX = set_colors(df_healthy, XKCD_COLORS_IDX,
                      DATA_FIGURE_COLORS, XKCD_COLORS)

    fig = plt.figure(figsize=(35,15))
    squeeze = 2
    gs = fig.add_gridspec(9,40 * squeeze, hspace = 0.75)

    axqpcr1 = fig.add_subplot(gs[2 : 4, 1*squeeze:14*squeeze])
    axrel1 = fig.add_subplot(gs[4:8,1*squeeze:14*squeeze])
    axpert1 = fig.add_subplot(gs[2:8,1*squeeze:14*squeeze], facecolor='none')
    axinoculum1 = fig.add_subplot(gs[4:8,0])

    #print("color index", XKCD_COLORS_IDX)
    max_qpcr_value_healthy, XKCD_COLORS_IDX, order_, order_inoc = plot_rel_and_qpcr(subjset_healthy,
    subjset_inoc = subjset_inoc, df = df_healthy, dset_type = "healthy",
    axqpcr = axqpcr1, axrel = axrel1, axpert=axpert1, axinoculum = axinoculum1,
    make_ylabels=True, color_taxa_dict = DATA_FIGURE_COLORS,
    color_index = XKCD_COLORS_IDX, color_set = XKCD_COLORS, figlabelinoculum
    = 'D', figlabelqpcr='B', figlabelrel='E', make_legend=False,
    taxaname_map = taxa_map_healthy, inoc_order = None)
    #print("color index", XKCD_COLORS_IDX)

    axqpcr2 = fig.add_subplot(gs[2 : 4,17 * squeeze : 30 * squeeze])
    axrel2 = fig.add_subplot(gs[4 : 8, 17 * squeeze : 30 * squeeze])
    axpert2 = fig.add_subplot(gs[2 : 8, 17 * squeeze : 30 * squeeze],
               facecolor='none')
    axinoculum2 = fig.add_subplot(gs[4 : 8, 16 * squeeze])

    max_qpcr_value_uc, XKCD_COLORS_IDX, order_, order_inoc = plot_rel_and_qpcr(subjset_uc,
    subjset_inoc = subjset_inoc, df = df_uc, dset_type = "uc", axqpcr = axqpcr2,
    axrel = axrel2, axpert=axpert2, axinoculum = axinoculum2, make_ylabels = True,
    color_taxa_dict = DATA_FIGURE_COLORS, color_index = XKCD_COLORS_IDX,
    color_set = XKCD_COLORS, figlabelinoculum ='F', figlabelqpcr = 'C',
    figlabelrel='G', make_legend = False, taxaname_map = taxa_map_uc,
    labels_order = order_, inoc_order = order_inoc)

    #make the qpcr value (absolute value) axis labels consistent in both axes
    max_qpcr_value = np.max([max_qpcr_value_healthy, max_qpcr_value_uc])
    axqpcr1.set_ylim(bottom = 1e9, top = max_qpcr_value * (1.25))
    axqpcr2.set_ylim(bottom = 1e9, top = max_qpcr_value * (1.25))
    axqpcr1.set_yscale('log')
    axqpcr2.set_yscale('log')

    axqpcr1.set_yticks([1e10, 1e11])
    for tick in axqpcr1.yaxis.get_major_ticks():
        tick.label.set_fontsize(24)

    axqpcr2.set_yticks([1e10, 1e11])
    for tick in axqpcr2.yaxis.get_major_ticks():
        tick.label.set_fontsize(24)

    ax_experiment =  fig.add_subplot(gs[0, 5 * squeeze : 25 * squeeze],
                      facecolor='none')
    add_expt_figure(ax_experiment, subjset_healthy, figlabel = 'A')

    #make legend
    axlegend = fig.add_subplot(gs[2 : 8, 33 * squeeze: 39 * squeeze],
                facecolor='none')
    plot_legend(axlegend = axlegend, level = TAXLEVEL, cutoff = CUTOFF_FRAC_ABUNDANCE,
    color_taxa_dict = DATA_FIGURE_COLORS)

    axpert1.spines['top'].set_visible(False)
    axpert1.spines['bottom'].set_visible(False)
    axpert1.spines['left'].set_visible(False)
    axpert1.spines['right'].set_visible(False)
    axpert1.xaxis.set_major_locator(plt.NullLocator())
    axpert1.xaxis.set_minor_locator(plt.NullLocator())
    axpert1.yaxis.set_major_locator(plt.NullLocator())
    axpert1.yaxis.set_minor_locator(plt.NullLocator())

    axpert2.spines['top'].set_visible(False)
    axpert2.spines['bottom'].set_visible(False)
    axpert2.spines['left'].set_visible(False)
    axpert2.spines['right'].set_visible(False)
    axpert2.xaxis.set_major_locator(plt.NullLocator())
    axpert2.xaxis.set_minor_locator(plt.NullLocator())
    axpert2.yaxis.set_major_locator(plt.NullLocator())
    axpert2.yaxis.set_minor_locator(plt.NullLocator())


    fig.subplots_adjust(wspace = 0.58, left = 0.05, right = 0.95, top =  0.925,
    bottom = .075, hspace = 0.8)

    plt.savefig(SAVEPATH + "main_figure2_.png", dpi = 100, bbox_inches = "tight")
    plt.savefig(SAVEPATH + "main_figure2_.pdf", dpi = 400, bbox_inches = "tight")

def plot_legend(axlegend, level, cutoff, color_taxa_dict):
    """plots the legend"""

    taxidx = TAXLEVEL_REV_IDX[TAXLEVEL]
    upper_tax = TAXLEVEL_INTS[taxidx+1]
    lower_tax = TAXLEVEL_INTS[taxidx]

    #print(level, upper_tax, lower_tax)
    labels = list(color_taxa_dict.keys())
    labels.sort()
    print(labels)

    #others appears last in the legend
    last_label = None
    for label in labels:
        print(label)
        if len(label.split(' ')) > 2:
            last_label = label
            break
    if last_label is not None:
        labels.remove(last_label)
        labels.append(last_label)

    ims = []
    for label in labels:
        im, = axlegend.bar([0],[0], color= color_taxa_dict[label], label=label)
        ims.append(im)

    extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none',
    linewidth=0)
    legend_handle = [extra]
    legend_handle = legend_handle + ims
    extra_col = [extra]*(len(ims)+1)
    legend_handle = legend_handle + extra_col + extra_col

    empty_label = ''
    legend_labels = [empty_label]* (len(ims)+1) + ['$\\bf{' +
    upper_tax.capitalize() + '}$']
    for label in labels[:-1]:
        l1,_ = label.split(' ')
        if l1 == 'nan':
            l1 = 'Uncultured Clone'
        legend_labels = legend_labels + [l1.capitalize()]
    legend_labels = legend_labels + ['Other < {}%'.format(CUTOFF_FRAC_ABUNDANCE*100)]
    legend_labels = legend_labels + ['$\\bf{' + lower_tax.capitalize() + '}$']
    for label in labels[:-1]:
        _,l2 = label.split(' ')
        print(label, l2)
        if l2 == 'NA':
            l2 = ''
        legend_labels = legend_labels + [l2.capitalize()]

    axlegend.legend(legend_handle, legend_labels, ncol = 3, loc='upper center',
        fontsize = 22, columnspacing =0.01)

    axlegend.spines['top'].set_visible(False)
    axlegend.spines['bottom'].set_visible(False)
    axlegend.spines['left'].set_visible(False)
    axlegend.spines['right'].set_visible(False)
    axlegend.xaxis.set_major_locator(plt.NullLocator())
    axlegend.xaxis.set_minor_locator(plt.NullLocator())
    axlegend.yaxis.set_major_locator(plt.NullLocator())
    axlegend.yaxis.set_minor_locator(plt.NullLocator())


main()
