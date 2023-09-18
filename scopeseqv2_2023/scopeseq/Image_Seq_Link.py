import os
import fnmatch
import pandas as pd
import numpy as np

from scopeseq.Cell_Bead_Link import CellBeadLink
from scopeseq.method import patch_index_replace
from scopeseq.utils import pickle_dump, pickle_load


class ImageSeqLink:
    def __init__(self, project_folder):
        self._project_folder = project_folder
        self._seq_folder = project_folder + 'seq/'
        self._analysis_folder = project_folder + 'analysis/'
        self._obc_cell = None
        self._cell = None
        self._cell_id = None
        self._seq_num = None
        self.seqobc_cell = None
        self._linked_num = None

    def set_image(self, image_fn):
        """

        :param image_fn: file name - CellBeadLink object file
        :return:
        """
        # check params
        print('image_fn: ' + image_fn)
        if not os.path.exists(image_fn):
            raise ValueError(f'image_fn does not exist')

        image = pickle_load(image_fn)

        if not isinstance(image, CellBeadLink):
            raise TypeError(f'bead should be CellBeadLink type.')

        self._obc_cell = image.obc_cell
        self._cell = image.cell

    def set_sequencing(self, cell_id=None, use='index'):
        """

        :param cell_id: list, np.ndarray - cell index
                        file name - cell index
        :param use: category - 'index', use index column as cell index
                               'columns', use columns names as cell index, 'gid', 'gene' are removed.
        :return:
        """
        if isinstance(cell_id, (list, np.ndarray)):
            self._cell_id = list(cell_id)

        elif isinstance(cell_id, str):
            # check params
            if not os.path.exists(cell_id):
                raise ValueError(f'cell_id file does not exist.')
            if use not in ['index', 'columns']:
                raise ValueError(f'use should be either "index" or "columns".')

            if use == 'index':
                ids = pd.read_csv(cell_id, usecols=[0], header=None, sep='\t')
                ids = list(ids[0])
                # remove
                if ids[0] != ids[0]:
                    ids.remove(ids[0])
                self._cell_id = ids

            if use == 'columns':
                f = open(cell_id, 'r')
                ids = f.read().split('\t')
                f.close()
                # remove
                for x in ['gid', 'gene', '\n']:
                    if x in id:
                        ids.remove(x)
                self._cell_id = ids

        else:
            raise ValueError(f'cell_id should be list, ndarray, or str.')

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

    def seq_image_link(self, replace):
        """

        :param replace: dict - {image patch: sequencing i7 index}
        :return:
        """
        print("linking sequencing barcode to optical barcode ...")

        # map optical barcode to sequencing i7 index
        seq_index = patch_index_replace(self._obc_cell['obc'], replace)

        # get unique optical barcode
        values, index, counts = np.unique(seq_index, return_index=True, return_counts=True)
        unique_index = sorted(index[counts == 1])
        obc_cell = self._obc_cell.loc[unique_index]
        unique_seq_index = seq_index[unique_index]
        obc_cell.index = unique_seq_index

        # get matched sequencing barcode
        matched_mask = np.isin(self._cell_id, unique_seq_index)
        matched_id = np.array(self._cell_id)[matched_mask]

        self.seqobc_cell = obc_cell.loc[matched_id]
        self._seq_num = len(self._cell_id)
        self._linked_num = len(matched_id)

        print(f'get {self._seq_num} sequenced cells.'
              f'linked {self._linked_num} cells.')

    def generate_sc_measurements(self, seq_hist=None, seq_matrix=None):
        self._check_params()

        # get linked cell barcode
        linked_ids = list(self.seqobc_cell.index)

        # image
        linked_cellindex = self.seqobc_cell['cell_index']
        sc_image_matrix = self._cell.loc[linked_cellindex]
        sc_image_matrix.index = linked_ids
        sc_image_matrix['cell_index'] = linked_cellindex
        sc_image_matrix.to_csv(self._analysis_folder + 'sc_image.matrix.linked.txt', sep='\t')

        # seq
        if seq_hist is not None:
            # check params
            print('seq_hist: ' + seq_hist)
            if not os.path.exists(seq_hist):
                raise ValueError(f'seq_hist does not exist.')

            print('seq_matrix: ' + seq_matrix)
            if seq_matrix is None:
                raise ValueError(f'seq_matrix cannot be empty.')
            if not os.path.exists(seq_matrix):
                raise ValueError(f'seq_matrix does not exist.')

            # hist
            hist = pd.read_csv(seq_hist, index_col=0, sep='\t', header=None)
            sc_seq_hist = hist.loc[linked_ids]
            sc_seq_hist.to_csv(self._analysis_folder + 'sc_seq.hist.linked.txt', sep='\t', header=None)

            # matrix
            matrix = pd.read_csv(seq_matrix, sep='\t', header=None)
            matrix.columns = ['gid', 'gene'] + hist.index.to_list()
            sc_seq_matrix = matrix[['gid', 'gene'] + linked_ids]
            sc_seq_matrix.to_csv(self._analysis_folder + 'sc_seq.matrix.linked.txt', sep='\t', header=None, index=False)

    def merge(self, obj1, obj2):
        # get two image_seq_link objects
        link1 = pickle_load(obj1)
        link2 = pickle_load(obj2)
        seqobc = pd.concat([link1.seqobc_cell, link2.seqobc_cell])

        # get unique cell-bead pair
        cellindex = seqobc['cell_index']
        values_cell, index_cell = np.unique(cellindex, return_index=True)
        values_bead, index_bead = np.unique(seqobc.index, return_index=True)
        index = list(set.intersection(set(index_bead), set(index_cell)))
        seqobc = seqobc.iloc[index, :]

        self.seqobc_cell = seqobc
        self._linked_num = len(index)

        # get info
        print(f'get {self._linked_num} linked cells.')
