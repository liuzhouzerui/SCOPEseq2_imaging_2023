import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import joblib
import anndata
import loompy
from scopeseq.clustering import cptt, run_ploidy
from scopeseq.plot import plt_map
from sklearn.decomposition import PCA
from scipy.stats import zscore


class Cluster:
    def __init__(self, project_folder):
        self._project_folder = project_folder
        self._analysis_folder = project_folder + 'analysis/'
        self._anndata = None
        self._ploidy_chr = None
        self._ploidy_zscore = None
        self._malignant_mask = None
        self._linked_mask = None
        self._image = None
        self._pcs = None

    def set_cluster(self, mtx, pg, umap):
        self._check_params()

        print('mtx: ' + mtx)
        if not os.path.exists(mtx):
            raise ValueError(f'mtx does not exist')

        print('pg: ' + pg)
        if not os.path.exists(pg):
            raise ValueError(f'pg does not exist')

        print('umap: ' + umap)
        if not os.path.exists(umap):
            raise ValueError(f'umap does not exist')

        # read in count matrix
        self._anndata = anndata.read_loom(mtx, obs_names='', var_names='')
        self._anndata.obs['index'] = self._anndata.obs['batch'].astype('str') + '-' + \
            self._anndata.obs['sample'].astype('str') + '-' + \
            self._anndata.obs['CellID'].astype('str')
        # read in phenograph cluster
        self._anndata.obs['pgs'] = pd.read_csv(pg, sep='\t', header=None)[0].values
        # read in umap coordinates
        self._anndata.obsm['umap_emb'] = pd.read_csv(umap, sep='\t', header=None).to_numpy()

    def plt_sample(self):
        pdf_outfile = self._analysis_folder + 'umap_sample.pdf'
        with PdfPages(pdf_outfile) as pdf:
            # color umap by sample
            colors = ['red', 'blue', 'green', 'magenta', 'brown', 'cyan', 'black', 'orange', 'grey', 'darkgreen',
                      'yellow', 'tan', 'seagreen', 'fuchsia', 'gold', 'olive']
            sample = self._anndata.obs['batch'].astype('str') + '-' + \
                self._anndata.obs['sample'].astype('str')
            plt_map(pdf, x=self._anndata.obsm['umap_emb'], name='Sample', color=colors,
                    group=sample, t='category', s=2)

    def link_ploidy(self, ploidy_file, sample, batch):
        # read in ploidy
        ploidy = pd.read_csv(ploidy_file, sep='\t', index_col=0, header=None)
        ploidy.index = str(batch) + '-' + sample + '-' + ploidy.index.values

        # initialize
        if 'ploidy' not in self._anndata.obs.columns:
            self._anndata.obs['ploidy'] = None

        # assign ploidy
        index = self._anndata.obs['index'][np.isin(self._anndata.obs['index'], ploidy.index)]
        index_origin = self._anndata.obs.index
        self._anndata.obs.index = self._anndata.obs['index']
        self._anndata.obs.loc[index, 'ploidy'] = ploidy.loc[index][4]
        self._anndata.obs.index = index_origin

    def _plt_ploidy(self):
        pdf_outfile = self._analysis_folder + 'umap_ploidy.pdf'
        with PdfPages(pdf_outfile) as pdf:
            # color umap by ploidy
            mask = self._anndata.obs['ploidy'].notna()
            plt_map(pdf, x=self._anndata.obsm['umap_emb'][mask, :], name='ploidy',
                    color=self._anndata.obs['ploidy'][mask], cmap='coolwarm', t='float', s=2)

    def ploidy_pca(self, gmt_file):
        self._ploidy_chr, self._ploidy_zscore = run_ploidy(gmt_file, gids=self._anndata.var['Accession'],
                                                           matrix=self._anndata.X.T, chr_use=[6, 9])
        pca = PCA(n_components=3)
        pcs = pca.fit_transform(zscore(self._ploidy_chr.T, axis=1))
        self._pcs = pd.DataFrame(pcs, columns=['PC1', 'PC2', 'PC3'])

    def assign_cell_type(self, cell_type_ref):
        self._anndata.obs['cell_type'] = None
        for cell_type in cell_type_ref.keys():
            self._anndata.obs['cell_type'][np.isin(self._anndata.obs['pgs'], cell_type_ref[cell_type])] = \
                cell_type

        self._plt_cell_type(cell_type_ref)

    def _plt_cell_type(self, cell_type_ref):
        pdf_outfile = self._analysis_folder + 'umap_celltype.pdf'
        with PdfPages(pdf_outfile) as pdf:
            # color umap by cell type
            plt_map(pdf, x=self._anndata.obsm['umap_emb'], name='Cell Type',
                    color=mpl.cm.tab10.colors, group=self._anndata.obs['cell_type'], order=cell_type_ref.keys(),
                    t='category', s=3)

    def _check_params(self):
        # folders
        print('project_folder: ' + self._project_folder)
        if not os.path.exists(self._project_folder):
            raise ValueError(f'project_folder does not exist')

        print('analysis_folder: ' + self._analysis_folder)
        if not os.path.exists(self._analysis_folder):
            raise ValueError(f'analysis_folder does not exist')

    def link_image(self, image_file, sample, batch=None):
        """

        :param image_file: file name - linked single cell image matrix, output from the ImageSeqLink object
        :param sample: str - sample name
        :param batch: str - batch name
        :return:
        """
        image = pd.read_csv(image_file, sep='\t', index_col=0)
        image.index = batch + '-' + sample + '-' + image.index.values

        # initialize
        if self._linked_mask is None:
            self._linked_mask = np.zeros(self._anndata.n_obs, dtype=bool)
        if self._image is None:
            self._image = pd.DataFrame(index=self._anndata.obs['index'], columns=image.columns)

        index = self._anndata.obs['index'][np.isin(self._anndata.obs['index'], image.index)]
        self._linked_mask = np.logical_or(np.isin(self._anndata.obs['index'], image.index), self._linked_mask)

        self._image.loc[index, image.columns] = image.loc[index]

    def _plt_image(self, channels):
        pdf_outfile = self._analysis_folder + 'umap_image.pdf'
        with PdfPages(pdf_outfile) as pdf:
            # color umap by image
            for i in channels:
                mask = self._image[i+'_signal'].notna().values
                plt_map(pdf, x=self._anndata.obsm['umap_emb'][mask, :], name=i,
                        color=self._image[i+'_signal'].values[mask], cmap='Reds', t='float', s=2)

    def to_anndata(self):
        # anndata, cells x genes
        # export to loom based hdf5
        write_obsm_varm = False
        if self._anndata.obsm:
            write_obsm_varm = True
        self._anndata.obs.to_csv(self._project_folder + 'annotation.txt', sep='\t', index=False)
        self._anndata.write_loom(self._project_folder + 'count_matrix.loom', write_obsm_varm=write_obsm_varm)
        return self._anndata

    def plt_markergene(self, marker):
        matrix_norm = cptt(self._anndata.X.T)
        pdf_outfile = self._analysis_folder + 'umap_gene.pdf'
        with PdfPages(pdf_outfile) as pdf:
            for i in marker:
                if i in self._anndata.var['Gene'].values:
                    plt_map(pdf, x=self._anndata.obsm['umap_emb'], name=i,
                            color=matrix_norm[(self._anndata.var['Gene'] == i).values, :].tolist()[0],
                            cmap='Reds', t='float', s=0.5)
                else:
                    print('marker:'+i+' is not in gene list')
