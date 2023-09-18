import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

from numpy.random import choice
from random import random
from scipy.optimize import leastsq
from scipy.stats import ranksums, mannwhitneyu, pearsonr
from statsmodels.sandbox.stats.multicomp import multipletests

# from rpy2.robjects.packages import importr
# import rpy2.robjects as ro
# from rpy2.rinterface import RRuntimeWarning
import warnings

from scopeseq.utils import load_matrix, write_matrix


"""
fit distribution
"""


def cutoff(x):
    chist, cedge = np.histogram(np.log2(x[x > 0]), bins=50)
    md = np.median(cedge)
    mn1 = cedge[np.argmax(chist[cedge[0:-1] < md])]
    mn2 = cedge[np.argmax(chist[cedge[0:-1] > md]) + len(cedge[cedge <= md])]
    # the shortest bin between two mean as thresh
    thresh = cedge[np.argmin(chist[(cedge[0:-1] > mn1) & (cedge[0:-1] < mn2)]) + len(cedge[cedge <= mn1])]
    return 2 ** thresh


# fit single gaussian
def gaussian(data, params):
    c1, mn1, st1 = params
    return c1*np.exp(-0.5*(data-mn1)*(data-mn1)/(st1*st1))


def gaussian_res(params, vals, data):
    return vals - gaussian(data, params)


def fit_gaussian(x, name, plt_output, tail='lower', alpha=0.01,  x2=None, sample=None, var=1, bins=50):
    # initial
    if sample is None:
        sample = ['sample', 'control']

    # fit
    chist, cedge = np.histogram(np.log2(x[(~np.isnan(x)) & (x > 1)]), bins=bins)
    params = []
    if sum(params) == 0:
        peak1 = np.max(chist)
        mn1 = cedge[np.argmax(chist)]
        params = [peak1, mn1, var]
    print(params)
    fit = leastsq(gaussian_res, params, args=(chist, cedge[0:len(cedge) - 1]))
    print(fit[0])

    if alpha == 0.01:
        a = 2.57
    elif alpha == 0.05:
        a = 1.96
    else:
        raise ValueError(f'alpha should be either "0.01" or "0.05".')

    if tail == 'lower':
        thresh = fit[0][1] - a * fit[0][2]  # a=0.05, use 1.96; a=0.01 use 2.57
    elif tail == 'upper':
        thresh = fit[0][1] + a * fit[0][2]
    else:
        raise ValueError(f'tail should be either "lower" or "upper".')

    print(2 ** thresh)
    # plot
    if x2 is None:
        plt_hist_threshold(fit, chist, cedge, thresh, name, plt_output, model='gaussian')
    else:
        plt_hist_threshold_2(x2, x, sample, fit, chist, cedge, thresh, name, plt_output, model='gaussian')
    return 2 ** thresh, fit[0]


# fit double gaussian
def double_gaussian(data, params):
    c1, mn1, st1, c2, mn2, st2 = params
    return c1*np.exp(-0.5*(data-mn1)*(data-mn1)/(st1*st1)) + c2*np.exp(-0.5*(data-mn2)*(data-mn2)/(st2*st2))


def double_gaussian_res(params, vals, data):
    return vals - double_gaussian(data, params)


def fit_double_gaussian(x, name, plt_output, tail='lower', alpha=0.01, x2=None, sample=None, var=1, bins=50, pos=True):
    # initial
    if sample is None:
        sample = ['sample', 'control']
    if pos is True:
        chist, cedge = np.histogram(np.log2(x[(~np.isnan(x)) & (x > 1)]), bins=bins)
    else:
        chist, cedge = np.histogram(np.log2(x[(~np.isnan(x))]), bins=bins)
    params = []
    if sum(params) == 0:
        md = np.median(cedge)
        peak1 = np.max(chist[cedge[0:-1] < md])
        peak2 = np.max(chist[cedge[0:-1] > md])
        mn1 = cedge[np.argmax(chist[cedge[0:-1] < md])]
        mn2 = cedge[np.argmax(chist[cedge[0:-1] > md]) + len(cedge[cedge <= md])]
        params = [peak1, mn1, var, peak2, mn2, var]
    print(params)
    fit = leastsq(double_gaussian_res, params, args=(chist, cedge[0:len(cedge) - 1]))
    print(fit[0])

    if alpha == 0.01:
        a = 2.57
    elif alpha == 0.05:
        a = 1.96
    else:
        raise ValueError(f'alpha should be either "0.01" or "0.05".')

    if tail == 'lower':
        thresh = fit[0][1] + a * fit[0][2]
    elif tail == 'upper':
        thresh = fit[0][4] - a * fit[0][5]
    else:
        raise ValueError(f'tail should be either "lower" or "upper".')

    print(2 ** thresh)
    # plot
    if x2 is None:
        plt_hist_threshold(fit, chist, cedge, thresh, name, plt_output, model='double_gaussian')
    else:
        plt_hist_threshold_2(x2, x, sample, fit, chist, cedge, thresh, name, plt_output, model='double_gaussian')
    return 2 ** thresh, fit[0]


def plt_hist_threshold(fit, chist, cedge, thresh, name, plt_output, model='gaussian'):
    # plot
    plt.bar(cedge[0:len(cedge) - 1], chist, color='k', width=(cedge[1] - cedge[0]) * 0.8)
    if model == 'gaussian':
        plt.plot(cedge[0:len(cedge)], gaussian(cedge[0:len(cedge)], fit[0]), 'r')
    if model == 'double_gaussian':
        plt.plot(cedge[0:len(cedge)], double_gaussian(cedge[0:len(cedge)], fit[0]), 'r')
    plt.vlines(x=thresh, ymin=0, ymax=chist.max())
    plt.xlabel('log2(' + name + ')')
    plt.ylabel('Number of Cells')
    plt.savefig(plt_output)
    plt.close()


def plt_hist_threshold_2(x1, x2, sample, fit, chist, cedge, thresh, name, plt_output, model='gaussian', color=None):
    if color is None:
        color = ['red', 'gray']

    plt.hist(np.log2(x1[(~np.isnan(x1)) & (x1 > 1)]), bins=50, color=color[0], alpha=0.5)
    plt.hist(np.log2(x2[(~np.isnan(x2)) & (x2 > 1)]), bins=50, color=color[1], alpha=0.5)
    plt.legend(sample)
    if model == 'gaussian':
        plt.plot(cedge[0:len(cedge)], gaussian(cedge[0:len(cedge)], fit[0]), 'r')
    if model == 'double_gaussian':
        plt.plot(cedge[0:len(cedge)], double_gaussian(cedge[0:len(cedge)], fit[0]), 'r')
    plt.vlines(x=thresh, ymin=0, ymax=chist.max())
    plt.xlabel('log2(' + name + ')')
    plt.ylabel('Number of Cells')
    plt.savefig(plt_output)
    plt.close()


"""
Mann-Whitney U test
"""


# def scran_normalize(matrixfile, samplefile, normfile):
#     Rcmd1 = 'library("scran");'
#     Rcmd2 = 'df <- read.table("%(matrixfile)s",sep="\t");' % vars()
#     Rcmd3 = 'data <- as.matrix(df[,3:ncol(df)]);' % vars()
#     Rcmd4 = 'clusters <- as.matrix(read.table("%(samplefile)s"));' % vars()
#     Rcmd5 = 'sfs <- as.vector(computeSumFactors(data,clusters=clusters));'
#     Rcmd6 = 'write.table(sfs,"%(normfile)s",sep="\t",row.names=FALSE,col.names=FALSE);' % vars()
#
#     warnings.filterwarnings("ignore", category=RRuntimeWarning)
#     ro.r(Rcmd1)
#     ro.r(Rcmd2)
#     ro.r(Rcmd3)
#     ro.r(Rcmd4)
#     ro.r(Rcmd5)
#     ro.r(Rcmd6)
#
#     return 0


def mwdiffex(filename, gids, genes, matrix1, matrix2):
    stats = []
    pvals = []
    lfcs = []
    mns1 = np.mean(matrix1,axis=1)
    mns2 = np.mean(matrix2,axis=1)
    i = 0
    for vec1, vec2, mn1, mn2 in zip(matrix1, matrix2, mns1, mns2):
        if mn1 * mn2 == 0.0:
            if mn1 == 0.0 and mn2 == 0.0:
                u = float('nan')
                p = -1.0
                lfcs.append(float('nan'))
            elif mn1 == 0.0:
                u, p = mannwhitneyu(vec1, vec2, alternative='two-sided')
                lfcs.append(float('inf'))
            elif mn2 == 0.0:
                u, p = mannwhitneyu(vec1, vec2, alternative='two-sided')
                lfcs.append(-1.0 * float('inf'))
        else:
            u, p = mannwhitneyu(vec1, vec2, alternative='two-sided')
            lfcs.append(np.log2(mn2/mn1))
        stats.append(u)
        pvals.append(p)
        i += 1
    stats = np.array(stats)
    pvals = np.array(pvals)
    lfcs = np.array(lfcs)
    qvals1 = multipletests(pvals[pvals != -1.0], method='fdr_bh')[1]
    qvals = []
    j = 0
    for p in pvals:
        if p != -1.0:
            qvals.append(qvals1[j])
            j += 1
        else:
            qvals.append(float('nan'))
    qvals = np.array(qvals)
    ind = np.argsort(qvals)
    stats = stats[ind]
    pvals = pvals[ind]
    lfcs = lfcs[ind]
    gids = gids[ind]
    genes = genes[ind]
    qvals = qvals[ind]
    mns1 = mns1[ind]
    mns2 = mns2[ind]

    with open(filename, 'w') as g:
        g.write('GID\tgene\tlog2(fold-change)\tstat\tsample1 mean\tsample2 mean\tpval\tqval\n')
        for gid, gene, lfc, stat, mn1, mn2, pval, qval in zip(gids, genes, lfcs, stats, mns1, mns2, pvals, qvals):
            g.write('%(gid)s\t%(gene)s\t%(lfc)f\t%(stat)f\t%(mn1)f\t%(mn2)f\t%(pval)e\t%(qval)e\n' % vars())

    return 0


def downsample_vector(vector, frac):
    molecules = []
    N = float(sum(vector))
    M = len(vector)
    for i in range(M):
        molecules.extend([i for m in range(vector[i])])
    dsmolecules = choice(molecules, int(frac*N), replace=False)
    dsvector = [np.count_nonzero(dsmolecules == i) for i in range(M)]
    return dsvector


def downsample_merged_matrix(matrix, samples):
    smasks = [samples==s for s in set(samples)]
    mncts = [np.mean(matrix[:, smask]) for smask in smasks]
    mincts = np.min(mncts)
    fracs = [float(mincts)/float(mnct) for mnct in mncts]
    newmatrix = []
    for s,vector in zip(samples, matrix.T):
        if fracs[s] < 1.0:
            newmatrix.append(downsample_vector(vector, fracs[s]))
        else:
            newmatrix.append(vector)
    return np.transpose(np.array(newmatrix))


def subsample_cells_matrix(matrix, Ncells):
    N = matrix.shape[1]
    rmask = np.zeros(N,dtype=int)
    rmask[:Ncells] = 1
    np.random.shuffle(rmask)
    rmask = rmask.astype(bool)
    return matrix[:, rmask]


def mwdiffex_prep(sample, count_matrix, output_dir, prefix, whitelist=None, blacklist=None):
    r = str(np.random.randint(0, 10000000))
    dsmatrix_outfile = output_dir + '/' + prefix + '.' + r + '.ds.matrix.txt'
    dssample_outfile = output_dir + '/' + prefix + '.' + r + '.ds.samples.txt'

    gids, genes, matrix = load_matrix(count_matrix, 0)
    if whitelist:
        wl = set([line.split()[0].split('.')[0] for line in open(whitelist)])
        wlmask = np.array([True if gid.split('.')[0] in wl else False for gid in gids])
        gids = gids[wlmask]
        genes = genes[wlmask]
        matrix = matrix[wlmask]
    if blacklist:
        bl = set([line.split()[0].split('.')[0] for line in open(blacklist)])
        blmask = np.array([True if gid.split('.')[0] in bl else False for gid in gids])
        gids = gids[~blmask]
        genes = genes[~blmask]
        matrix = matrix[~blmask]

    unmask = np.array([True if s == 0 else False for s in sample])
    unmatrix = matrix[:, unmask]
    N_u = np.sum(unmask)
    trmask = np.array([True if s == 1 else False for s in sample])
    trmatrix = matrix[:, trmask]
    N_t = np.sum(trmask)
    if N_t > N_u:
        trmatrix = subsample_cells_matrix(trmatrix, N_u)
        mmatrix = np.concatenate((unmatrix, trmatrix), axis=1)
        msamples = np.array([0 if i < N_u else 1 for i in range(N_u + N_u)])
    elif N_u > N_t:
        new_unmatrix = subsample_cells_matrix(unmatrix, N_t)
        mmatrix = np.concatenate((new_unmatrix, trmatrix), axis=1)
        msamples = np.array([0 if i < N_t else 1 for i in range(N_t + N_t)])
    else:
        mmatrix = np.concatenate((unmatrix, trmatrix), axis=1)
        msamples = np.array([0 if i < N_t else 1 for i in range(N_t + N_t)])
    if len(msamples) > 40:
        mmatrix = downsample_merged_matrix(mmatrix, msamples)
        write_matrix(dsmatrix_outfile, gids, genes, mmatrix)
        np.savetxt(dssample_outfile, msamples, fmt='%d')
    return r, gids, genes, mmatrix, msamples


def mwdiffex_run(output_dir, prefix, r, gids, genes, mmatrix, msamples, scran=False):
    dsmatrix_outfile = output_dir + '/' + prefix + '.' + r + '.ds.matrix.txt'
    dssample_outfile = output_dir + '/' + prefix + '.' + r + '.ds.samples.txt'
    norm_outfile = output_dir + '/' + prefix + '.' + r + '.ds.norm.txt'
    if scran:
        scran_normalize(dsmatrix_outfile, dssample_outfile, norm_outfile)
    norm = np.loadtxt(norm_outfile)
    minnz = np.min(norm[norm > 0.0])
    norm = np.array([n if n > minnz else minnz for n in norm])
    mmatrix = mmatrix / norm
    msamples = msamples.astype(bool)
    new_unmatrix = mmatrix[:, ~msamples]
    trmatrix = mmatrix[:, msamples]
    wc_outfile = output_dir + '/' + prefix + '.mwdiffex.txt'
    mwdiffex(wc_outfile, gids, genes, new_unmatrix, trmatrix)


