import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


def plt_hist(x, name, output):
    plt.hist(np.log2(x[(~np.isnan(x)) & (x > 1)]), bins=50, color='k')
    plt.xlabel('log2(' + name + ')')
    plt.ylabel('Number of Cells')
    output.savefig()
    plt.close()


def plt_two_hist(x1, x2, sample, name, output, color=None):
    if color is None:
        color = ['red', 'gray']

    plt.hist(np.log2(x1[(~np.isnan(x1)) & (x1 > 1)]), bins=50, color=color[0], alpha=0.5)
    plt.hist(np.log2(x2[(~np.isnan(x2)) & (x2 > 1)]), bins=50, color=color[1], alpha=0.5)
    plt.legend(sample)
    plt.xlabel('log2(' + name + ')')
    plt.ylabel('Number of Cells')
    output.savefig()
    plt.close()


def plt_map(pdf, x, name, color, group=None, order=None, s=5, cmap='coolwarm', t='float', **kwargs):
    if t == 'float':
        mask = [False if np.isnan(i) else True for i in color]
        na_mask = [not i for i in mask]
        fig, ax = plt.subplots()
        ax.scatter(x[:, 0][na_mask], x[:, 1][na_mask], c='gray', alpha=0.1, s=s, lw=0, edgecolor='none')
        im = ax.scatter(x[:, 0][mask], x[:, 1][mask], c=np.array(color)[mask], s=s, lw=0, edgecolor='none', cmap=cmap,
                        **kwargs)
        ax.set_xlabel('Dim 1')
        ax.set_ylabel('Dim 2')
        ax.set_title('Colored by ' + name)
        plt.colorbar(im, ax=ax)
        pdf.savefig(bbox_inches='tight')
        plt.close()

    elif t == 'category':
        if order is None:
            cat = sorted(set(group))
        else:
            cat = order

        fig, ax = plt.subplots()
        for i, c in enumerate(cat):
            x_sub = x[group == c, :]
            im = ax.scatter(x_sub[:, 0], x_sub[:, 1], color=color[i], s=s, lw=0, edgecolor='none', label=c, **kwargs)
        ax.legend(loc='center left', bbox_to_anchor=(1.2, 0.5))
        ax.set_xlabel('Dim 1')
        ax.set_ylabel('Dim 2')
        ax.set_title('Colored by ' + name)
        plt.colorbar(im, ax=ax)
        pdf.savefig(bbox_inches='tight')
        plt.close()
