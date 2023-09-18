import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

from scipy.stats import gmean
from scopeseq.plot import plt_map


"""
# hashtag
htos = ['Hashtag1','Hashtag2','Hashtag3']
samples = ['DMSO','Etoposide','Panobinostat']
hashtag = hto(project_folder,prefix,htos,samples)
"""


def CLRnorm(vec):
    return np.log1p(vec/gmean(vec[vec>0]))


def hto(project_folder, prefix, htos, samples, umap_file=None):
    # read in hashtag counts
    hashtag = pd.read_csv(project_folder+'seq/'+prefix+'_HTO.demultiplexed.txt', sep='\t', index_col=0)
    hto_sample_map = pd.Series(samples, index=htos)
    hto_sample_map.loc['Doublet']='Doublet'
    hashtag['sample'] = hashtag['Classification'].apply(lambda x: hto_sample_map[x])

    # umap
    # read in UMAP coordinates
    if umap_file:
        umap = pd.read_csv(umap_file, sep='\t', header=None)
        umap.index = hashtag.index

    # Normalize and classify each barcode for plotting
    hashtag_norm = hashtag[htos].apply(lambda x: CLRnorm(x), axis=0)
    hashtag_norm['sample'] = hashtag['sample']
    hashtag_frac = hashtag[htos].divide(hashtag[htos].sum(axis=1), axis=0)

    pdf_outfile = project_folder + 'analysis/hashtag.pdf'
    with PdfPages(pdf_outfile) as pdf:
        data_plt = hashtag['Classification'].value_counts()
        data_plt = data_plt.reindex(htos+['Doublet'])
        plt.bar(data_plt.index, data_plt.values)
        plt.ylabel('Number of cells')
        pdf.savefig(bbox_inches='tight')
        plt.close()

        # grouping and plotting cells by negative, singlet, doublet
        for i in htos:
            plt.scatter(hashtag_norm[i], hashtag_frac[i], s=1, c=(hashtag['Classification']!=i)*1, cmap='nipy_spectral')
            plt.title(i)
            plt.xlabel('Normalized Counts')
            plt.ylabel('Fraction of Raw Counts')
            pdf.savefig(bbox_inches='tight')
            plt.close()

        for i in htos:
            plt.hist(hashtag_norm[i], bins=30)
            plt.title(i)
            plt.xlabel('Normalized Counts')
            plt.ylabel('Cell counts')
            pdf.savefig()
            plt.close()

    if umap_file:
        pdf_outfile = project_folder + 'analysis/umap_hto.pdf'
        with PdfPages(pdf_outfile) as pdf:
            for i in htos:
                plt.scatter(umap[0], umap[1], c=hashtag_norm[i], s=5, cmap='Reds')
                plt.xlabel('UMAP 1')
                plt.ylabel('UMAP 2')
                plt.colorbar()
                plt.title('colored by '+i+' normalized counts')
                pdf.savefig(bbox_inches='tight')
                plt.close()

            plt_map(pdf, x=umap.values, name='colored by sample',
                        color=mpl.cm.tab20.colors, group=hashtag['sample'].values, order=samples+['Doublet'],
                        t='category', s=5)

    # order data
    hashtag_norm_order = pd.DataFrame()
    for i in range(len(samples)):
        data = hashtag_norm[hashtag['sample']==samples[i]]
        data = data.sort_values(hashtag_norm.columns[i], ascending=False)
        hashtag_norm_order = pd.concat([hashtag_norm_order, data])
    data = hashtag_norm[hashtag['sample']=='Doublet']
    hashtag_norm_order = pd.concat([hashtag_norm_order, data])

    lut = dict(zip(samples+['Doublet'], sns.color_palette("tab20")[0:(len(htos)+1)]))
    row_colors = hashtag_norm_order['sample'].map(lut)

    pdf_outfile = project_folder + 'analysis/hashtag.heatmap.pdf'
    with PdfPages(pdf_outfile) as pdf:
        sns.clustermap(hashtag_norm_order[htos], row_cluster=False, col_cluster=False,
                       row_colors = row_colors,
                       cmap='Reds', yticklabels=False)
        pdf.savefig(bbox_inches='tight')
        plt.close()

    # export
    hashtag.to_csv(project_folder+prefix+'_HTO.demultiplexed.txt', sep='\t')
    return hashtag