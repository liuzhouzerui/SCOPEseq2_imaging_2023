import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import anndata
import loompy
from scopeseq.utils import load_matrix


class CountMatrix:
    def __init__(self, project_folder):
        self._project_folder = project_folder
        self._analysis_folder = project_folder + 'analysis/'
        self._matrix_file = None
        self._hist_file = None
        self._matrix = None
        self._n_gene = None
        self._n_cell = None
        self._gene_summary = None
        self._cell_summary = None
        self._sample = None

    def set_project_folder(self, project_folder):
        self._project_folder = project_folder
        self._analysis_folder = project_folder + 'analysis/'

    def _set_matrix(self, matrix_file):
        print('set matrix...')

        # get count matrix, gene information
        self._matrix_file = matrix_file
        gids, genes, self._matrix = load_matrix(matrix_file, subsample=0)

        self._n_gene = self._matrix.shape[0]
        self._n_cell = self._matrix.shape[1]

        self._gene_summary = pd.DataFrame({'Accession': gids, 'Gene': genes})
        self._cell_summary = pd.DataFrame(None, index=range(self._n_cell))

    def _set_hist(self, hist_file):
        print('set hist...')

        if self._cell_summary is None:
            raise RuntimeError(f'count matrix does not exist. set matrix first.')

        # get cell information
        self._hist_file = hist_file
        hist = pd.read_csv(hist_file, sep='\t', header=None, names=['CellID', 'molecules', 'genes'])
        self._cell_summary = pd.concat([self._cell_summary, hist], axis=1)

    def _check_mt(self):
        print('set mt...')

        # mitochondria
        mt = self._matrix[self._gene_summary['Gene'].str.contains('^MT-'), :].sum(axis=0)
        self._cell_summary['mt'] = mt
        self._cell_summary['mt_percentage'] = mt/self._cell_summary['molecules']

    def _set_sample(self, seq_index):
        print('set sample...')

        if not isinstance(seq_index, dict):
            raise TypeError(f'seq_index should be dict type.')

        if 'CellID' not in self._cell_summary.columns:
            raise ValueError(f'cell_summary does not have CellID column. set hist first.')

        self._sample = list(seq_index.keys())

        self._cell_summary['sample'] = None

        for sample, ids in seq_index.items():
            if seq_index[sample]:
                cell_id = self._cell_summary['CellID'].apply(lambda x: int(x.split('_')[0]))
                self._cell_summary.loc[np.isin(cell_id, ids), 'sample'] = sample
            else:
                self._cell_summary['sample'] = sample

    def _set_batch(self, seq_index):
        print('set batch...')

        if not isinstance(seq_index, dict):
            raise TypeError(f'seq_index should be dict type.')

        self._batch = list(seq_index.keys())

        self._cell_summary['batch'] = None

        for batch, ids in seq_index.items():
            if seq_index[batch]:
                cell_id = self._cell_summary['CellID'].apply(lambda x: int(x.split('_')[0]))
                self._cell_summary.loc[np.isin(cell_id, ids), 'batch'] = batch
            else:
                self._cell_summary['batch'] = batch

    def set_sequencing(self, matrix_file, hist_file, sample_index, batch_index):
        """

        :param matrix_file: count matrix file name - no header, gene x cell. first column: gene id, second column: gene name
        :param hist_file: sequencing history file name - no header. cell x ['CellID', 'molecules', 'genes']
        :param sample_index: dict - {sample name: sequencing i7 index list}
        :param batch_index: dict - {batch name: sequencing i7 index list}
        :return:
        """
        self._set_matrix(matrix_file)
        self._set_hist(hist_file)
        self._check_mt()
        self._set_sample(sample_index)
        self._set_batch(batch_index)

    def seq_summary(self):
        self._check_params()

        self._cell_summary.to_csv(self._analysis_folder + 'seq_summary.txt', sep='\t', index=False)

        pdf_outfile = self._analysis_folder + 'seq_summary.pdf'
        with PdfPages(pdf_outfile) as pdf:
            for col in self._cell_summary.columns:
                if pd.api.types.is_numeric_dtype(self._cell_summary[col]):
                    sns.violinplot(x='sample', y=col, data=self._cell_summary)
                    plt.xlabel('Sample')
                    plt.ylabel(col)
                    pdf.savefig()
                    plt.close()

        # get info
        cell_num = self._cell_summary.shape[0]
        print(f'get {cell_num} sequenced cells.')

    def _check_params(self):
        # folders
        print('project_folder: ' + self._project_folder)
        if not os.path.exists(self._project_folder):
            raise ValueError(f'project_folder does not exist')

        print('analysis_folder: ' + self._analysis_folder)
        if not os.path.exists(self._analysis_folder):
            raise ValueError(f'analysis_folder does not exist')

    def to_anndata(self):
        # anndata, cells x genes
        data = anndata.AnnData(X=self._matrix.T, obs=self._cell_summary, var=self._gene_summary)

        # export to loom based hdf5
        data.write_loom(self._project_folder + 'count_matrix.loom')

        return data
