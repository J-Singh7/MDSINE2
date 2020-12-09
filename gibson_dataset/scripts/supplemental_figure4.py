"""
python supplemental_figure4.py -file1 "../../processed_data/gibson_healthy_agg_taxa.pkl" \
       -file2 "../../processed_data/gibson_uc_agg_taxa.pkl" \
       -file3 "../../processed_data/gibson_inoculum_agg_taxa.pkl"
"""

import mdsine2 as md2

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

import seaborn as sns

import numpy as np
import pickle

from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests, fdrcorrection
import skbio
import skbio.diversity
import skbio.stats.distance
import skbio.diversity
#import pylab as pl
import argparse


HEALTHY_SUBJECTS = ['2','3','4','5']
UNHEALTHY_SUBJECTS = ['6','7','8','9','10']
PERTURBATION_COLOR = 'orange'


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

def add_perturbation_label(ax, perturbations, subj, times, textcolor='black',
     textsize=None, alpha=0.25, label=True):
    """
       adds perturbation labels to the axes, ax
    """
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
        pert_locs.append((times.index(perturbation.ends[subj]) + times.index(perturbation.starts[subj])) / 2)
        name = perturbation.name
        if name is None:
            name = 'pert{}'.format(pidx)
        pert_names.append(name)

    if label:
        # Set the names on the top x-axis
        ax2 = ax.twiny()

        # # Set the visibility of the twin axis to see through
        # ax2.spines['top'].set_visible(False)
        # ax2.spines['bottom'].set_visible(False)
        # ax2.spines['left'].set_visible(False)
        # ax2.spines['right'].set_visible(False)
        # ax2.xaxis.set_major_locator(plt.NullLocator())
        # ax2.xaxis.set_minor_locator(plt.NullLocator())
        # ax2.yaxis.set_major_locator(plt.NullLocator())
        # ax2.yaxis.set_minor_locator(plt.NullLocator())

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


def _set_type_of_point_bc(t, subj, subjset):
    '''Sets the type of marker that is used in beta diversity figure
    each perturbation, post perturbation, and colonization gets a different number.

    Here we know there are three perturbations
    '''

    if t < 21.5:
        return 0
    elif t >= 21.8 and t <= 28.5:
        return 1
    elif t > 28.5 and t < 35.5:
        return 2
    elif t >= 35.5 and t <= 42.5:
        return 3
    elif t > 42.5 and t < 50.5:
        return 4
    elif t >= 50.5 and t <= 57.5:
        return 5
    else:
        return 6

def unbias_var_estimate(vals):
    '''
       @parameters
       -----------------------------------------------------------------------
       vals : (list[float])

       @returns
       (float) Unbiased variance estimate of the values in vals
    '''

    vals = np.asarray(vals)
    mean = np.mean(vals)
    a = np.sum((vals - mean)**2)
    return a / (len(vals)-1)

def compute_alpha_diversity_time(subjset):
    """
       computes the alpha diversity over time for each subject

       @parameters
       -----------------------------------------------------------------------
       subjset : ([pylab.subject])

       @returns
       -----------------------------------------------------------------------
       (dict (float) -> [float])
    """

    div_values = {}
    for subj in subjset:
        for t in subj.times:
            val = md2.diversity.alpha.normalized_entropy(subj.reads[t])
            if t not in div_values:
                div_values[t] = []
            div_values[t].append(val)

    return div_values

def compute_mean_std(alpha_values_d, time_li):

    """
       @parameters
       -----------------------------------------------------------------------
       alpha_values_d : (dict (float) -> [float]) alpha diversity values for
                        each subject at each time
       time_li : ([time]) times at which samples are collected

       @returns
       -----------------------------------------------------------------------
       (lists of floats) : Mean and Standard deviation in alpha_values_d
    """

    mean_li = []
    std_li = []
    for i, t in enumerate(time_li):
        mean_li.append(np.mean(alpha_values_d[t]))
        std_li.append(np.sqrt(unbias_var_estimate(alpha_values_d[t])))

    return np.asarray(mean_li), np.asarray(std_li)

def compute_p_alpha_diversity(alpha_healthy_d, alpha_uc_d, time_li):
    """
       computes the p-values using Mann Whitney test

       @parameters
       -----------------------------------------------------------------------
       alpha_healthy_d : (dict (float) -> [float])
       alpha_uc_d : (dict (float) -> [float])
       time_li : ([float])

       @returns
       -----------------------------------------------------------------------
       results of the Mann-Whitney test (tuple consisting of hypothesis test
       results and the adjusted p-values)
    """

    p_vals = []
    for t in time_li:
        test_res = mannwhitneyu(alpha_healthy_d[t], alpha_uc_d[t],
                   use_continuity = False)
        p_vals.append(test_res[1])

    hypothesis_test_result = multipletests(p_vals, alpha = 0.05,
    method = "fdr_bh", is_sorted = False)

    return hypothesis_test_result

def classify_p_values(time_array, p_details):
    """
       classify the p-values as significant / not significant and return the
       the two types of values separately

       @parameters
       -----------------------------------------------------------------------
       time_array : ([float]) list of times
       p_details : (tuple) results of the Mann-Whitney test

       @returns
       -----------------------------------------------------------------------
       list of significant times, list of significant p-values, lists of
       non-significant times, list of non-significant p-values
    """

    time_1 = []
    time_2 = []
    p_1 = []
    p_2 = []
    for i in range(len(time_array)):
        if p_details[0][i] == True:
            time_1.append(i)
            p_1.append(p_details[1][i])
        else:
            time_2.append(i)
            p_2.append(p_details[1][i])

    return time_1, p_1, time_2, p_2

def alpha_diversity_mean_std(subjset_healthy, subjset_uc, subjset_innoc, name,
    ax = None, axlegend = None, figlabel = None, save = False):

    """
       Plots the mean alpha diversity values and the standard deviation

       @parameters
       -----------------------------------------------------------------------
       subjset_heallthy, subjset_uc : ([pylab.subject]) data pertaining to
                                       healthy / uc subjsets
       subjset_innoc : ([pylab.subject]) data at the time of inoculum
       name : (str) figure name
       ax, axlegend : (matplotlib.Axes)
       save : (bool) save figure or not
       figlabel : (str) figure label

       @returns
       -----------------------------------------------------------------------
       None
    """

    print("Plotting Alpha Diversity")
    #subjset_healthy = md2.Study.load('../gibson_output/datasets/gibson_healthy_agg.pkl')
    #subjset_uc = md2.Study.load('../gibson_output/datasets/gibson_uc_agg.pkl')
    #subjset_innoc = md2.Study.load('../gibson_output/datasets/gibson_inoculum_agg.pkl')

    if ax is None:
        fig = plt.figure(figsize = (12, 8))
        ax = fig.add_subplot(1, 1, 1)

    alpha_healthy = compute_alpha_diversity_time(subjset_healthy)
    alpha_uc = compute_alpha_diversity_time(subjset_uc)
    alpha_innoc_healthy = md2.diversity.alpha.normalized_entropy(
         subjset_innoc['Healthy'].reads[0])
    alpha_innoc_uc = md2.diversity.alpha.normalized_entropy(
             subjset_innoc['Ulcerative Colitis'].reads[0])
    subj_ = ""
    for subj in subjset_healthy:
        subj_ = subj
        break
    name_ = subj_.name

    times = np.sort(list(alpha_healthy.keys()))
    means_healthy, std_healthy = compute_mean_std(alpha_healthy, times)
    means_uc, std_uc = compute_mean_std(alpha_uc, times)

    test_results = compute_p_alpha_diversity(alpha_healthy, alpha_uc, times)
    sig_times, sig_p, non_sig_times, non_sig_p = classify_p_values(times,
     test_results)

    times_idxs = np.arange(len(times))
    #print('times:', times)

    times_idxs_healthy = times_idxs-(0.25/2)
    times_idxs_uc = times_idxs+(0.25/2)

    colors = sns.color_palette('muted')
    colors_healthy = colors[0]
    colors_uc = colors[1]

    #plot the error bar
    ax.errorbar(times_idxs_healthy, means_healthy, std_healthy,
        ecolor=colors_healthy, color=colors_healthy, capsize=3, fmt='none')
    ax.plot(times_idxs_healthy, means_healthy, marker='o', color=colors_healthy,
    linewidth=0, markersize=5)

    ax.errorbar(times_idxs_uc, means_uc, std_uc,
        ecolor=colors_uc, color=colors_uc, capsize=3, fmt='none')
    ax.plot(times_idxs_uc, means_uc, marker='o', color=colors_uc, linewidth=0,
        markersize=5)

    # Add in the Inoculum
    ax.plot([-1.5], alpha_innoc_healthy, marker='*', color=colors_healthy,
    markersize=10)
    ax.plot([-1.5], alpha_innoc_uc, marker='*', color=colors_uc, markersize=10)

    #plot markers inidicating whether p-values are significant/non-significant
    ax.plot(sig_times, [0.2] * len(sig_times), marker = "x", linewidth = 0,
    color = "red")
    ax.plot(non_sig_times, [0.2] * len(non_sig_times), marker = "o",
    linewidth = 0, color = "red")

    # Set the xticklabels
    locs = np.arange(0, len(times),step=10)
    ticklabels = times[locs]
    plt.xticks(ticks=locs, labels=ticklabels)

    # Plot the perturbation markers
    for perturbation in subjset_healthy.perturbations:
        start_t = perturbation.starts['2']
        end_t = perturbation.ends['2']
        print(start_t)
        print(end_t)

        startidx = np.searchsorted(times, start_t)
        endidx = np.searchsorted(times, end_t)
        perturbation.start = startidx - 0.5
        perturbation.end = endidx + 0.5
    times_li = list(times)
    add_perturbation_label(ax, subjset_healthy.perturbations,subj_, times_li,
    textsize=18, alpha=0)
    for perturbation in subjset_healthy.perturbations:
        ax.axvline(x = times_li.index(perturbation.starts[name_]), color='black',
        linestyle='--', lw=2)
        ax.axvline(x = times_li.index(perturbation.ends[name_]), color='black',
        linestyle='--', linewidth=2)

    # Make legend
    if axlegend is None:
        axlegend = fig.add_subplot(111, facecolor='none')
        axlegend.spines['top'].set_visible(False)
        axlegend.spines['bottom'].set_visible(False)
        axlegend.spines['left'].set_visible(False)
        axlegend.spines['right'].set_visible(False)
        axlegend.xaxis.set_major_locator(plt.NullLocator())
        axlegend.xaxis.set_minor_locator(plt.NullLocator())
        axlegend.yaxis.set_major_locator(plt.NullLocator())
        axlegend.yaxis.set_minor_locator(plt.NullLocator())

    handles = []
    l = mlines.Line2D([],[], color=colors_healthy,
        linestyle='-', label='Healthy', linewidth = 5)
    handles.append(l)
    l = mlines.Line2D([],[], color=colors_uc,
        linestyle='-', label='Ulcerative Colitis', linewidth = 5)
    handles.append(l)
    l = mlines.Line2D([], [], color='black', marker='*',
        label='Inoculum', linestyle='none', markersize=10)
    handles.append(l)
    lgnd2 = axlegend.legend(handles=handles, fontsize=17, title_fontsize = 18,
        title='$\\bf{Dataset}$', bbox_to_anchor=(1.04, 1.025),
        loc = "upper left")
    axlegend.add_artist(lgnd2)

    handles = []
    l = mlines.Line2D([],[], color="red", linestyle='none',
        marker='x', label='Significant')
    handles.append(l)
    l = mlines.Line2D([],[], color="red",linestyle='none',
        marker ='o', label='Not significant')
    handles.append(l)
    lgnd3 = axlegend.legend(handles = handles, title='$\\bf{P-values}$',
            loc='lower left', borderaxespad=0., bbox_to_anchor=(1.05, 0),
            title_fontsize=18, fontsize = 17)
    axlegend.add_artist(lgnd3)

    # Set the ticks to be bold
    ax.set_yticks([0.2, 0.4, 0.6, 0.8])
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(18)
        # tick.label.set_fontweight('bold')
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(18)
        # tick.label.set_fontweight('bold')

    # Set the labels
    ax.set_ylabel('nat', size=20, fontweight='bold')
    ax.set_xlabel('Time (days)', size=20, fontweight='bold')
    ax.set_title("A", fontsize = 25, fontweight = 'bold', loc = "left")

    plt.subplots_adjust(bottom=0.18)

    # caption = 'Mean and standard deviation of normalized entropy measure within each consortium.'
    # axlegend.text(0.5, -0.2, caption, horizontalalignment='center', fontsize=17)
    if save:
        plt.savefig(name + 'alpha_diversity_mean_std.pdf', bbox_inches = "tight")
        plt.savefig(name + 'alpha_diversity_mean_std.png', bbox_inches = "tight")
    # plt.show()
    plt.close()

def reformat_data(subjset):
    """
        format the data to make it compatible for computing Bray-Curtis
        dissimilarity

        @parameters
        ----------------------------------------------------------------------
        subjset : (pl.base.Subject)

        @returns
        ----------------------------------------------------------------------
        reformated data (numpy), sample labels ([str]),
        sample labels id ((dict (float) -> (str, float)))
    """

    data = None
    labels = []
    labels_float = {}
    for subj in subjset:
        print(subj.name)
        for t in subj.times:
            ts = str(float(t)).replace('.5', 'PM').replace('.0', 'AM')
            labels.append('{}-{}'.format(subj.name, ts))
            labels_float[labels[-1]] = (subj, t)
        d = subj.matrix()["raw"].T
        if data is None:
            data = d
        else:
            data = np.vstack((data, d))

    print()
    return data, labels, labels_float

def combine_dicts(dict1, dict2):

    """
       combines two dictionaries into one single dictionary
    """

    dict_ = {}
    for keys1 in dict1:
        dict_[keys1] = dict1[keys1]

    for keys2 in dict2:
        dict_[keys2] = dict2[keys2]

    return dict_

def beta_diversity_figure(subjset_healthy, subjset_uc, subjset_innoc, name = None,
    axleft = None, axright = None, axcenter = None, figlabel = None,
    save = False):

    """
       Plots the first two dimensions of the beta diversity

       @parameters
       -----------------------------------------------------------------------
       subjset_heallthy, subjset_uc : ([pylab.subject]) data pertaining to
                                       healthy / uc subjsets
       subjset_innoc : ([pylab.subject]) data at the time of inoculum
       name : (str) figure name
       ax, axleft, axright, axcenter : (matplotlib.Axes)
       save : (bool) save figure or not
       figlabel : (str) figure label

       @returns
       -----------------------------------------------------------------------
       None
    """

    print("Running Beta Diversity")
    heal_data, heal_labels, heal_labels_float = reformat_data(subjset_healthy)
    uc_data, uc_labels, uc_labels_float = reformat_data(subjset_uc)

    print(uc_data.shape, heal_data.shape)
    data = np.vstack((heal_data, uc_data))
    labels = heal_labels + uc_labels
    labels_float = combine_dicts(heal_labels_float, uc_labels_float)

    for subj in subjset_innoc:
        labels.append('inoculum {}'.format(subj.name))
        #print(subj.name)
        labels_float[labels[-1]] = 1000
        m = subj.matrix()['raw'].T
        #print(m.shape)
        data = np.vstack((data, m))

    print(data.shape, len(labels), len(labels_float))
    bc_dm = skbio.diversity.beta_diversity("braycurtis", data, labels)
    bc_pcoa = skbio.stats.ordination.pcoa(bc_dm)

    data = bc_pcoa.samples.to_numpy()

    if axleft is None:
        fig = plt.figure(figsize=(16,8))
        axright = fig.add_subplot(121)
        axleft = fig.add_subplot(122)

    colors = sns.color_palette('muted')
    colorshealthy = colors[0]
    colorsuc = colors[1]
    colorinoculumhealthy = colors[0]
    colorinoculumuc = colors[1]
    x_healthy = None
    y_healthy = None
    x_uc = None
    y_uc = None

    xs = []
    ys = []
    cs = []

    subj_colors = {}
    subj_perts = {}

    for subj in subjset_healthy:
        subj_colors[subj.name] = colorshealthy
        subj_perts[subj.name] = [[[],[]] for i in range(7)]

    for subj in subjset_uc:
        subj_colors[subj.name] = colorsuc
        subj_perts[subj.name] = [[[],[]] for i in range(7)]

    for row in range(data.shape[0]):
        if 'inoculum' in labels[row]:
            if 'Healthy' in labels[row]:
                print("healthy innoc")
                x_healthy = -data[row,0]
                y_healthy = data[row,1]
            else:
                print("uc innoc")
                x_uc = -data[row,0]
                y_uc = data[row,1]
        else:
            subj,t = labels_float[labels[row]]
            mi = _set_type_of_point_bc(t, subj, subjset_healthy)

            subj_perts[subj.name][mi][0].append(-data[row,0])
            subj_perts[subj.name][mi][1].append(data[row,1])

    PERT_MARKERS = ['+', 'd', 'o', 's', 'v', 'x', 'X']
    # PERT_MARKERS = ['+', 'd', 'v', 's', 'v', 'x', 'v']
    INOCULUM_MARKER = '*'
    for ix, ax in enumerate([axright, axleft]):
        print('ix', ix)

        # Plot the points
        for subj in subjset_healthy:
            for mi in range(len(subj_perts[subj.name])):
                xs = subj_perts[subj.name][mi][0]
                ys = subj_perts[subj.name][mi][1]

                ax.plot(xs, ys, PERT_MARKERS[mi], color=subj_colors[subj.name],
                    markersize=6, alpha=0.75)
        for subj in subjset_uc:
            for mi in range(len(subj_perts[subj.name])):
                xs = subj_perts[subj.name][mi][0]
                ys = subj_perts[subj.name][mi][1]

                ax.plot(xs, ys, PERT_MARKERS[mi], color=subj_colors[subj.name],
                    markersize=6, alpha=0.75)

        # Plot the inoculum
        print(x_healthy, y_healthy, x_uc, y_uc)
        ax.plot([x_healthy], [y_healthy], INOCULUM_MARKER, color=colorinoculumhealthy,
            markersize=15, alpha=0.75)
        ax.plot([x_uc], [y_uc], INOCULUM_MARKER, color=colorinoculumuc,
            markersize=15, alpha=0.75)
    axleft.set_xlim(left=-.35, right=0)
    axleft.set_ylim(bottom=-0.18, top=-0.18+0.35)

    if axcenter is None:
        axcenter = fig.add_subplot(111, facecolor='none')
    axcenter.spines['top'].set_visible(False)
    axcenter.spines['bottom'].set_visible(False)
    axcenter.spines['left'].set_visible(False)
    axcenter.spines['right'].set_visible(False)
    axcenter.xaxis.set_major_locator(plt.NullLocator())
    axcenter.xaxis.set_minor_locator(plt.NullLocator())
    axcenter.yaxis.set_major_locator(plt.NullLocator())
    axcenter.yaxis.set_minor_locator(plt.NullLocator())

    if figlabel is not None:
        axcenter.text(x = -0.1, y= 1.1,
            s=figlabel, fontsize=30, fontweight='bold',
            transform=axcenter.transAxes)

    # Make the legend for the right side
    # Ulcerative colitis
    handles = []
    l = mlines.Line2D([],[], color=colorshealthy,
        linestyle='-', label='Healthy', linewidth=5)
    handles.append(l)
    l = mlines.Line2D([],[], color=colorsuc,
        linestyle='-', label='Ulcerative Colitis', linewidth=5)
    handles.append(l)
    lgnd2 = axcenter.legend(handles=handles, title='$\\bf{Dataset}$',
        bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,
        fontsize=17, title_fontsize=18)
    axcenter.add_artist(lgnd2)

    # Perturbations
    handles = []
    l = mlines.Line2D([],[],
        marker=INOCULUM_MARKER,
        markeredgecolor='black',
        markerfacecolor='none',
        linestyle='none',
        markersize=15,
        label='Inoculum')
    handles.append(l)
    l = mlines.Line2D([],[],
        marker=PERT_MARKERS[0],
        markeredgecolor='black',
        markerfacecolor='none',
        linestyle='none',
        markersize=15,
        label='Colonization')
    handles.append(l)
    #for perturbation in subjset_healthy.perturbations:
    #    print(perturbation.name)
    pert_names = []
    for p in subjset_healthy.perturbations:
        pert_names.append(p.name)
        pert_names.append(p.name)
    print(pert_names)
    for pidx in range(1,7):
        pert_name = pert_names[pidx - 1]

        # If pidx-1 is odd, then it is post perturbation
        if (pidx-1)%2 == 1:
            pert_name = 'Post-' + pert_name

        l = mlines.Line2D([],[],
            marker=PERT_MARKERS[pidx],
            markeredgecolor='black',
            markerfacecolor='none',
            linestyle='none',
            markersize=15,
            label=pert_name)
        handles.append(l)
    lgnd3 = axcenter.legend(handles=handles, title='$\\bf{Markers}$',
        bbox_to_anchor=(1.05, 0.0), loc='lower left', borderaxespad=0.,
        fontsize=17, title_fontsize=18)
    axcenter.add_artist(lgnd3)

    axcenter.set_xlabel('PC1: {:.3f}'.format(bc_pcoa.proportion_explained[0]),
        fontsize=20, fontweight='bold')
    axcenter.xaxis.set_label_coords(0.5,-0.08)

    axright.set_ylabel('PC2: {:.3f}'.format(bc_pcoa.proportion_explained[1]),
        fontsize=20, fontweight='bold')

    # mark_inset(parent_axes=axright, inset_axes=axleft, loc1a=2, loc1b=1,
    #     loc2a=3, loc2b=4, fc='none', ec='crimson')
    axleft.spines['top'].set_color('gray')
    axleft.spines['bottom'].set_color('gray')
    axleft.spines['left'].set_color('gray')
    axleft.spines['right'].set_color('gray')
    axleft.tick_params(axis='both', color='gray')

    rect =[patches.Rectangle(xy=(-0.35, -0.18), width=0.35, height=.35)]

    pc = PatchCollection(rect, facecolor='none', alpha=0.8, edgecolor='gray')
    axright.add_collection(pc)

    axleft.set_xticks([-0.35, -0.25, -0.15, -0.05])
    axleft.set_yticks([-0.15, -0.05, 0.05, 0.15])
    for tick in axleft.yaxis.get_major_ticks():
        tick.label.set_fontsize(18)
        # tick.label.set_fontweight('bold')
    for tick in axleft.xaxis.get_major_ticks():
        tick.label.set_fontsize(18)
        # tick.label.set_fontweight('bold')

    axright.set_yticks([-0.4, -0.2, 0, 0.2, 0.4])
    for tick in axright.yaxis.get_major_ticks():
        tick.label.set_fontsize(18)
        # tick.label.set_fontweight('bold')
    for tick in axright.xaxis.get_major_ticks():
        tick.label.set_fontsize(18)
        # tick.label.set_fontweight('bold')

    #fig.subplots_adjust(left=0.09, right=0.775, top=0.86)

    axright.set_title('B', fontsize=25, fontweight='bold', loc = "left")
    axleft.set_title('C', fontsize=25, fontweight='bold',loc = "left")

    # # Label the left and right side
    # axcenter.text(x=0., y=1.01, s='A', fontsize=25, fontweight='bold')
    # axcenter.text(x=0.545, y=1.01, s='B', fontsize=25, fontweight='bold')

    # bc_pcoa.plot()
    if save:
        plt.savefig(name + 'pcoa_braycurtis_w_zoom.pdf', bbox_inches = "tight")
        plt.savefig(name + 'pcoa_braycurtis_w_zoom.png', bbox_inches = "tight")
    plt.close()


def get_reads(subjset, type_):

    """
       gets read counts, groups and label names associated with each samples

       @Parameters
       -----------------------------------------------------------------------
       subjset : (pl.base.Subject)
       type_ : (str) healthy or UC cohort

       @returns
       -----------------------------------------------------------------------
       read at each samples (numpy), label names (list), cohort names (list)
    """

    colonization = -1

    reads = []
    groups = []
    labels = []
    i = 0
    for subjidx, subj in enumerate(subjset):
        for t in subj.times:
            if t < colonization:
                continue
            reads.append(subj.reads[t])
            cohort = type_
            if t < subjset.perturbations[0].start:
                area = 'post colonization'
            elif t >= subjset.perturbations[0].start and t <=subjset.perturbations[0].end:
                area = subjset.perturbations[0].name
            elif t > subjset.perturbations[0].end and t < subjset.perturbations[1].start:
                area = 'post {}'.format(subjset.perturbations[0].name)
            elif t >= subjset.perturbations[1].start and t <= subjset.perturbations[1].end:
                area = subjset.perturbations[1].name
            elif t > subjset.perturbations[1].end and t < subjset.perturbations[2].start:
                area = 'post {}'.format(subjset.perturbations[1].name)
            elif t >= subjset.perturbations[2].start and t <= subjset.perturbations[2].end:
                area = subjset.perturbations[2].name
            else:
                area = 'post {}'.format(subjset.perturbations[2].name)
            label = '{}-{}-{}'.format(cohort,area,i)
            print(label)
            i += 1
            group = '-'.join(label.split('-')[:-1])
            labels.append(label)
            groups.append(cohort)
            print(group)

    return np.asarray(reads), groups, labels

def permanova(subjset_healthy, subjset_uc):

    """
       runs the permanova test and prints the result

       @Parameters
       -----------------------------------------------------------------------
       subjset_healthy, subjset_uc : (pl.base.Subject)

       @returns
       -----------------------------------------------------------------------
       None

    """
    colonization = 5
    i = 0

    healthy_reads, healthy_group, healthy_labels = get_reads(subjset_healthy, "healthy")
    uc_reads, uc_group, uc_labels = get_reads(subjset_uc, "uc")
    print(healthy_reads.shape, len(healthy_group))
    print(uc_reads.shape, len(uc_group))
    reads = np.vstack((healthy_reads, uc_reads))
    #print(reads.shape)
    #print(reads[0])
    #print(healthy_reads[0])
    grouping = healthy_group + uc_group
    labels = healthy_labels + uc_labels
    print()
    bc_dm = skbio.diversity.beta_diversity(counts = reads, ids = labels , metric="braycurtis")
    print(bc_dm.shape, len(grouping))
    test_result = skbio.stats.distance.permanova(distance_matrix = bc_dm, grouping = grouping)
    print(test_result)

def diversity_plot(subjset_healthy, subjset_uc, subjset_innoc, name = None):

    """
       plots the alpha and beta diversity together
       @parameters
       -----------------------------------------------------------------------
       subjset_heallthy, subjset_uc : ([pylab.subject]) data pertaining to
                                       healthy / uc subjsets
       subjset_innoc : ([pylab.subject]) data at the time of inoculum
    """

    SAVEPATH = "output_figures/"
    fig = plt.figure(figsize = (15, 15))
    spec = GridSpec(nrows = 2, ncols = 2, hspace = 0.25)

    ax2 = fig.add_subplot(spec[1, 0 : 2])
    ax3 = fig.add_subplot(spec[1, 0])
    ax4 = fig.add_subplot(spec[1, 1])
    ax1 = fig.add_subplot(spec[0, 0 : 2])
    #ax5 = fig.add_subplot(spec[0, 0 : 2])

    beta_diversity_figure(subjset_healthy, subjset_uc, subjset_innoc,
    name = name, axleft = ax4, axright = ax3, axcenter = ax2)
    alpha_diversity_mean_std(subjset_healthy, subjset_uc, subjset_innoc,
    name = name, ax = ax1, axlegend = ax1)
    #alpha_diversity_mean_std()
    fig.savefig(SAVEPATH + "supplemental_figure4.pdf", bbox_inches = "tight")
    plt.close()

def loaddata(healthy):
    '''Load the subjectset data

    Parameters
    ----------
    healthy : bool, None
        If True, return the healthy subjects
        If False, return the non-healthy subjects
        If None, return all of the subjects

    Returns
    -------
    pl.base.SubjectSet
    '''
    subjset = pl.base.SubjectSet.load('../datasets/real_subjectset.pkl')
    # subjset.pop_times([0])
    if healthy is not None:
        if not pl.isbool(healthy):
            raise TypeError('`healthy` ({}) must be a bool'.format(healthy))
        if healthy:
            subjset.pop_subject(UNHEALTHY_SUBJECTS)
        else:
            subjset.pop_subject(HEALTHY_SUBJECTS)
    return subjset

def main():

    args = parse_args()
    subjset_healthy = md2.Study.load(args.healthy_pkl)
    subjset_uc = md2.Study.load(args.uc_pkl)
    subjset_inoc = md2.Study.load(args.inoc_pkl)


    #subjset_healthy = loaddata(True)
    #subjset_uc = loaddata(False)
    #subjset_innoc = pl.base.SubjectSet.load('../datasets/inoculum_subjectset.pkl')

    #to_delete = []
    #for asv in subjset_innoc.asvs:
    #    if asv.name not in subjset_healthy.asvs:
    #        to_delete.append(asv.name)
    #print("len:", len(to_delete))
    #subjset_innoc.pop_asvs(to_delete)

    #alpha_diversity_mean_std(subjset_healthy, subjset_uc, subjset_innoc, "new_", save = True)
    #beta_diversity_figure(subjset_healthy, subjset_uc, subjset_innoc, "new_", save= True)
    #permanova(subjset_healthy, subjset_uc)
    #beta_diversity_figure()
    diversity_plot(subjset_healthy, subjset_uc, subjset_inoc)

main()
