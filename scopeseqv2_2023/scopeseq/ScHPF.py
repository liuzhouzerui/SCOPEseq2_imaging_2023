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
from scopeseq.clustering import run_phenograph, run_umap, run_ploidy
from scopeseq.plot import plt_map


class ScHPF:
    def __init__(self, project_folder):
        self._project_folder = project_folder
        self._schpf_folder = project_folder + 'schpf/'
        self._analysis_folder = project_folder + 'analysis/'
        self._schpf = None
        self._schpf_anndata = None
        self._image = None
        self._factor_num = None
        self._pg_n = None
        self._pgcolors = None
        self._imputed_expression = None
        self._ploidy_chr = None
        self._ploidy_zscore = None
        self._malignant_mask = None
        self._linked_mask = None
        self._live_mask = None
        self._image_th = {}

    def set_schpf(self, mtx, genes, schpf_joblib):
        """
        initialize schpf model
        :param mtx: file name - filtered count matrix
        :param genes: file name - filtered gene list
        :param schpf_joblib: file name - schpf joblib file
        :return:
        """
        print('mtx: ' + mtx)
        if not os.path.exists(mtx):
            raise ValueError(f'mtx does not exist')

        print('genes: ' + genes)
        if not os.path.exists(genes):
            raise ValueError(f'genes does not exist')

        print('schpf_joblib: ' + schpf_joblib)
        if not os.path.exists(schpf_joblib):
            raise ValueError(f'schpf_joblib does not exist')

        # read in schpf filtered count matrix
        self._schpf_anndata = anndata.read_mtx(mtx)
        # read in schpf filtered gene list
        self._schpf_anndata.var = pd.read_csv(genes, usecols=[0, 1], sep='\t', header=None, names=['Accession', 'Gene'])
        # read in schpf model
        self._schpf = joblib.load(schpf_joblib)
        # compute cell score
        self._schpf_anndata.obsm['cell_score'] = self._schpf.cell_score()
        # compute gene score
        self._schpf_anndata.varm['gene_score'] = self._schpf.gene_score()
        # get number of factors in the schpf model
        self._factor_num = self._schpf_anndata.obsm['cell_score'].shape[1]

    def set_cell_annotation(self, annotation):
        """

        :param annotation: file name - cell annotation file. cell x features
        :return:
        """
        print('annotation: ' + annotation)
        if not os.path.exists(annotation):
            raise ValueError(f'annotation does not exist')
        # read in cell annotation
        obs = pd.read_csv(annotation, sep='\t')

        # check 'batch','sample','CellID' are in the anndata.obs column names
        if sum(np.isin(['batch', 'sample', 'CellID'], obs.columns)) != 3:
            raise ValueError(f'batch or sample or CellID is not in the anndata.obs column names')

        obs['index'] = obs['batch'].astype('str') + '-' + \
            obs['sample'].astype('str') + '-' + \
            obs['CellID'].astype('str')

        # add cell anntation to schpf filtered count matrix
        obs.index = self._schpf_anndata.obs.index
        self._schpf_anndata.obs = self._schpf_anndata.obs.join(obs)

    def plt_cell_score(self):
        if 'sample' not in self._schpf_anndata.obs.columns:
            raise ValueError(f'cannot find sample annotation. set sample through set_cell_annotation first.')

        if 'cell_score' not in self._schpf_anndata.obsm.keys():
            raise ValueError(f'cannot find cell_score observation. set schpf first.')

        sample = self._schpf_anndata.obs['batch'].astype('str') + '-' + self._schpf_anndata.obs['sample'].astype('str')

        pdf_outfile = self._analysis_folder + 'schpf_cell_score.pdf'
        with PdfPages(pdf_outfile) as pdf:
            for i in range(self._factor_num):
                sns.boxplot(x=sample, y=self._schpf_anndata.obsm['cell_score'][:, i])
                plt.title('factor ' + str(i))
                plt.xlabel('Sample')
                plt.ylabel('Cell Score')
                pdf.savefig()
                plt.close()

    def set_batch_factor(self, batch_factor=None):
        """

        :param batch_factor: int - schpf factor removed from downstream analysis
        :return:
        """
        if batch_factor is None:
            batch_factor = []

        self._schpf_anndata.obsm['cs_batch'] = np.delete(self._schpf_anndata.obsm['cell_score'], batch_factor, axis=1)
        self._schpf_anndata.varm['gs_batch'] = np.delete(self._schpf_anndata.varm['gene_score'], batch_factor, axis=1)

    def get_phenograph(self):
        self._schpf_anndata.obs['pgs'], self._pg_n, self._pgcolors = \
            run_phenograph(self._schpf_anndata.obsm['cs_batch'])

        # output
        self._schpf_anndata.obs['pgs'].to_csv(self._schpf_folder+'schpf.pgs.txt', header=False, index=False, sep='\t')

    def _plt_phenograph(self):
        if self._pg_n is not None:
            pdf_outfile = self._analysis_folder + 'schpf_umap_cluster.pdf'
            with PdfPages(pdf_outfile) as pdf:
                # color umap by cell cluster
                for i in range(self._pg_n):
                    mask = self._schpf_anndata.obs['pgs'] == i
                    non_mask = self._schpf_anndata.obs['pgs'] != i
                    plt.scatter(self._schpf_anndata.obsm['umap_emb'][non_mask, 0],
                                self._schpf_anndata.obsm['umap_emb'][non_mask, 1], c='gray', s=1)
                    plt.scatter(self._schpf_anndata.obsm['umap_emb'][mask, 0],
                                self._schpf_anndata.obsm['umap_emb'][mask, 1], c='red', s=1)
                    plt.title(str(i))
                    pdf.savefig(bbox_inches='tight')
                    plt.close()

    def get_umap(self):
        self._schpf_anndata.obsm['umap_emb'] = run_umap(self._schpf_anndata.obsm['cs_batch'])

        # output
        umap = pd.DataFrame(self._schpf_anndata.obsm['umap_emb'])
        umap.to_csv(self._schpf_folder + 'schpf.umap.txt', header=False, index=False, sep='\t')

    def plt_umap(self):
        pdf_outfile = self._analysis_folder + 'schpf_umap.pdf'
        with PdfPages(pdf_outfile) as pdf:
            # color umap by cluster
            plt_map(pdf, x=self._schpf_anndata.obsm['umap_emb'], name='Cluster', color=self._pgcolors,
                    group=self._schpf_anndata.obs['pgs'], t='category')

            # color umap by sample
            colors = ['red', 'blue', 'green', 'magenta', 'brown', 'cyan', 'black', 'orange', 'grey', 'darkgreen',
                      'yellow', 'tan', 'seagreen', 'fuchsia', 'gold', 'olive']
            sample = self._schpf_anndata.obs['batch'].astype('str') + '-' + \
                self._schpf_anndata.obs['sample'].astype('str')
            plt_map(pdf, x=self._schpf_anndata.obsm['umap_emb'], name='Sample', color=colors,
                    group=sample, t='category')

            # color umap by obs
            for col in self._schpf_anndata.obs.columns:
                if pd.api.types.is_numeric_dtype(self._schpf_anndata.obs[col]):
                    plt_map(pdf, x=self._schpf_anndata.obsm['umap_emb'], name=col,
                            color=np.log2(self._schpf_anndata.obs[col]).to_list(), cmap='coolwarm', t='float')

            # color umap by cell score
            for i in range(self._factor_num):
                plt_map(pdf, x=self._schpf_anndata.obsm['umap_emb'], name='Cell Score, Factor ' + str(i),
                        color=np.log2(self._schpf_anndata.obsm['cell_score'][:, i]).tolist(), cmap='coolwarm',
                        t='float')

        # color umap by cluster
        self._plt_phenograph()

    def impute_count_matrix(self):
        # imputed matrix
        self._imputed_expression = np.dot(self._schpf.beta.e_x, self._schpf.theta.e_x.T)

    def generate_ploidy(self, gmt_file, chr_use=None, amplify=None):
        # define ploidy chromosomes
        if chr_use is None:
            chr_use = [6, 9]
        if amplify is None:
            amplify = [True, False]

        # calculate ploidy score
        self._ploidy_chr, self._ploidy_zscore = run_ploidy(gmt_file, gids=self._schpf_anndata.var['Accession'],
                                                           matrix=self._imputed_expression, chr_use=chr_use, amplify=amplify)

    def get_malignant(self, mscore_thresh, frac_thresh=0.9):
        # malignant fraction in cluster
        mfracs = []
        for pg in set(self._schpf_anndata.obs['pgs']):
            if pg > -1:
                pmask = np.isin(self._schpf_anndata.obs['pgs'], pg)
                pscore = self._ploidy_zscore[pmask]
                mfracs.append(float(len(pscore[pscore > mscore_thresh])) / float(len(pscore)))

        mstatus = []
        for i, p in enumerate(self._schpf_anndata.obs['pgs']):
            if p > -1 and mfracs[p] > frac_thresh:
                mstatus.append(1)
            else:
                mstatus.append(0)
        mstatus = np.array(mstatus)

        self._malignant_mask = np.isin(mstatus, 1)

        ploidy_outfile = self._analysis_folder + 'malignant_ratio.imputed.txt'
        ploidy_out = np.column_stack((self._schpf_anndata.obs['CellID'], self._ploidy_chr[6, :], self._ploidy_chr[9, :],
                                      self._ploidy_zscore, mstatus))
        np.savetxt(ploidy_outfile, ploidy_out, delimiter='\t', fmt='%s\t%f\t%f\t%f\t%d')

        self._plt_ploidy()

    def _plt_ploidy(self):
        pdf_outfile = self._analysis_folder + 'schpf_umap_ploidy.pdf'
        with PdfPages(pdf_outfile) as pdf:
            # hist
            plt.hist(self._ploidy_zscore, bins=50)
            plt.xlabel('ploidy zscore')
            plt.ylabel('count')
            pdf.savefig()
            plt.close()

            # violinplot
            sns.violinplot(x=self._schpf_anndata.obs['pgs'], y=self._ploidy_zscore)
            plt.xlabel('phenograph cluster')
            plt.ylabel('ploidy zscore')
            pdf.savefig()
            plt.close()

            # color umap by ploidy zscore
            plt_map(pdf, x=self._schpf_anndata.obsm['umap_emb'], name='ploidy zscore',
                    color=self._ploidy_zscore, cmap='coolwarm', t='float')

            # color umap by ploidy
            plt_map(pdf, x=self._schpf_anndata.obsm['umap_emb'], name='ploidy',
                    color=self._malignant_mask*1, cmap='coolwarm', t='float')

    def assign_cell_type(self, cell_type_ref):
        """

        :param cell_type_ref: a directory. {'cell type': [clusters]}
        :return:
        """
        self._schpf_anndata.obs['cell_type'] = None
        for cell_type in cell_type_ref.keys():
            self._schpf_anndata.obs['cell_type'][np.isin(self._schpf_anndata.obs['pgs'], cell_type_ref[cell_type])] = \
                cell_type

        self._plt_cell_type(cell_type_ref)

    def _plt_cell_type(self, cell_type_ref):
        pdf_outfile = self._analysis_folder + 'schpf_umap_celltype.pdf'
        with PdfPages(pdf_outfile) as pdf:
            # color umap by ploidy zscore
            plt_map(pdf, x=self._schpf_anndata.obsm['umap_emb'], name='Cell Type',
                    color=mpl.cm.tab10.colors, group=self._schpf_anndata.obs['cell_type'], order=cell_type_ref.keys(),
                    t='category')

    def _check_params(self):
        # folders
        print('project_folder: ' + self._project_folder)
        if not os.path.exists(self._project_folder):
            raise ValueError(f'project_folder does not exist')

        print('seq_folder: ' + self._seq_folder)
        if not os.path.exists(self._seq_folder):
            raise ValueError(f'seq_folder does not exist')

        print('analysis_folder: ' + self._analysis_folder)
        if not os.path.exists(self._analysis_folder):
            raise ValueError(f'analysis_folder does not exist')

    def link_image(self, image_file, sample=None, batch=None):
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
            self._linked_mask = np.zeros(self._schpf_anndata.n_obs, dtype=bool)
        if self._image is None:
            self._image = pd.DataFrame(index=self._schpf_anndata.obs['index'], columns=image.columns)

        index = self._schpf_anndata.obs['index'][np.isin(self._schpf_anndata.obs['index'], image.index)]
        self._linked_mask = np.logical_or(np.isin(self._schpf_anndata.obs['index'], image.index), self._linked_mask)

        self._image.loc[index, image.columns] = image.loc[index]

    def image_cluster(self, sample, batch, thres):
        """

        :param sample: str - sample name
        :param batch: str - batch name
        :param thres: dict - {fluorescent channel: threshold}
        :return:
        """
        # initialize
        if 'image_cluster' not in self._schpf_anndata.obs.columns:
            self._schpf_anndata.obs['image_cluster'] = ''

        self._image_th[batch+'_'+sample] = thres

        c = 0
        for i in thres.keys():
            pos = (self._schpf_anndata.obs['batch'].astype('str') == batch) & \
                  (self._schpf_anndata.obs['sample'].astype('str') == sample) & \
                  (self._image[i + '_signal'].values > thres[i])
            neg = (self._schpf_anndata.obs['batch'].astype('str') == batch) & \
                  (self._schpf_anndata.obs['sample'].astype('str') == sample) & \
                  (self._image[i + '_signal'].values <= thres[i])
            if c > 0:
                self._schpf_anndata.obs.loc[pos, 'image_cluster'] = self._schpf_anndata.obs.loc[
                                                                        pos, 'image_cluster'] + '_pos'
                self._schpf_anndata.obs.loc[neg, 'image_cluster'] = self._schpf_anndata.obs.loc[
                                                                        neg, 'image_cluster'] + '_neg'
            else:
                self._schpf_anndata.obs.loc[pos, 'image_cluster'] = 'pos'
                self._schpf_anndata.obs.loc[neg, 'image_cluster'] = 'neg'
            c += 1

    def to_anndata(self):
        """
        write schpf filterd count matrix
        :return:
        """
        # anndata, cells x genes
        # export to loom based hdf5
        self._schpf_anndata.write_loom(self._project_folder + 'schpf.loom', write_obsm_varm=True)

        return self._schpf_anndata

    def marker_genes(self, factors=None, n=100):
        """
        get top n marker genes of all the selected factors
        :param factors: None, will use batch_factors removed gene score, gs_batch.
                        Or a list, indicating witch factors to use.
        :param n: the number of marker genes per factor
        :return:
        """
        # get gene score of selected schpf factors
        if factors is None:
            gene_score = self._schpf_anndata.varm['gs_batch']
        else:
            gene_score = self._schpf_anndata.varm['gene_score']
            gene_score = gene_score[:, factors]
        # get top n genes of all selected factors
        ranks = np.argsort(gene_score, axis=0)[::-1]
        ranks_top = ranks[0:n, :].T.flatten()
        marker_genes = self._schpf_anndata.var.iloc[ranks_top]
        # output
        marker_genes.to_csv(self._schpf_folder + 'marker_genes.txt', header=False, index=False, sep='\t')

        return marker_genes

    def to_gsea(self):
        """
        generate rank file of each schpf factor as input of gsea analysis
        :return:
        """
        gene_score = self._schpf_anndata.varm['gene_score']
        genes = self._schpf_anndata.var['Gene']
        for i in range(gene_score.shape[1]):
            gsea = pd.concat([genes, pd.DataFrame(gene_score[:, i])], axis=1)
            gsea.to_csv(self._schpf_folder + str(i) + '.rnk', sep='\t', index=False, header=False)
