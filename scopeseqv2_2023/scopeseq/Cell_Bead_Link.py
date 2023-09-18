import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from scopeseq.Bead_Intensity import BeadIntensity
from scopeseq.Cell_Intensity import CellIntensity
from scopeseq.method import register, rotate_position
from scopeseq.utils import pickle_dump, pickle_load


class CellBeadLink:
    def __init__(self, project_folder):
        self._project_folder = project_folder
        self._analysis_folder = project_folder + 'analysis/'
        self._bead_position = None
        self._bead_num = None
        self._patch = None
        self._n_patch_per_lane = None
        self._n_lane = None
        self._well_position = None
        self._well_num = None
        self._d_th = None
        self._obc = None
        self.cell = None
        self._landmark_fn = None
        self._bead_rotated_position = None
        self._well_bead_link = None
        self.obc_cell = None

    def set_project_folder(self, project_folder):
        self._project_folder = project_folder
        self._analysis_folder = project_folder + 'analysis/'

    def set_bead_cell(self, bead_fn, cell_fn, landmark_fn):
        """

        :param bead_fn: file name - BeadIntensity object file
        :param cell_fn: file name - CellIntensity object file
        :param landmark_fn: file name - register markers of cell and bead scan
        :return:
        """
        # check params
        print('bead_fn: ' + bead_fn)
        if not os.path.exists(bead_fn):
            raise ValueError(f'bead_fn does not exist.')

        print('cell_fn: ' + cell_fn)
        if not os.path.exists(cell_fn):
            raise ValueError(f'cell_fn does not exist.')

        print('landmark_fn: ' + landmark_fn)
        if not os.path.exists(landmark_fn):
            raise ValueError(f'landmark_fn does not exist.')

        # get bead and cell intensity
        bead = pickle_load(bead_fn)
        cell = pickle_load(cell_fn)

        if not isinstance(bead, BeadIntensity):
            raise TypeError(f'bead should be BeadIntensity type.')

        if not isinstance(cell, CellIntensity):
            raise TypeError(f'cell should be CellIntensity type.')

        # check if bead matches cell
        if (not bead.n_lane == cell.n_lane) | (not bead.total_lanes == cell.total_lanes):
            raise ValueError(f'BeadIntensity and CellIntensity object does not match.')

        # get information
        self._bead_position = bead.bead_position
        self._bead_num = bead.bead_num
        self._patch = bead.patch
        self._n_patch_per_lane = bead.n_patch_per_lane
        self._n_lane = bead.n_lane
        self._well_position = cell.well_position
        self._well_num = cell.well_num
        self._d_th = cell.d_th
        self._obc = bead.obc
        self.cell = cell.cell
        self._landmark_fn = landmark_fn

    def _check_params(self):
        # folders
        print('project_folder: ' + self._project_folder)
        if not os.path.exists(self._project_folder):
            raise ValueError(f'project_folder does not exist')

        print('analysis_folder: ' + self._analysis_folder)
        if not os.path.exists(self._analysis_folder):
            raise ValueError(f'analysis_folder does not exist')

    def _rotate_bead_position(self):
        print('mapping bead scan to cell scan...')

        # rotate bead position to match well position
        self._bead_rotated_position = pd.DataFrame(-1, index=range(self._bead_num), columns=['XM', 'YM'])

        # get landmark
        try:
            landmark_all = pd.read_csv(self._landmark_fn)
        except:
            landmark_all = pd.read_excel(self._landmark_fn, sep='\t')

        if landmark_all.shape[0] == 4:
            i = 0

            position = self._bead_position

            # get landmark
            landmark = landmark_all.iloc[(i * 4):((i + 1) * 4), :]
            # get device position
            target_start = np.array(landmark.loc[0 + (i * 4), ['XM', 'YM']])
            target_end = np.array(landmark.loc[1 + (i * 4), ['XM', 'YM']])
            initial_start = np.array(landmark.loc[2 + (i * 4), ['XM', 'YM']])
            initial_end = np.array(landmark.loc[3 + (i * 4), ['XM', 'YM']])

            # get rotated bead position
            self._bead_rotated_position[['XM', 'YM']] = rotate_position(position, target_start, target_end,
                                                                        initial_start, initial_end)

        elif landmark_all.shape[0] == 4*self._n_patch_per_lane:
            # for each patch in the current lane
            for i in range(self._n_patch_per_lane):
                # get patch number
                patch = i + self._n_patch_per_lane * self._n_lane

                # get bead position
                position = self._bead_position[self._patch == patch]

                # get landmark
                landmark = landmark_all.iloc[(i * 4):((i + 1) * 4), :]
                # get device position
                target_start = np.array(landmark.loc[0 + (i * 4), ['XM', 'YM']])
                target_end = np.array(landmark.loc[1 + (i * 4), ['XM', 'YM']])
                initial_start = np.array(landmark.loc[2 + (i * 4), ['XM', 'YM']])
                initial_end = np.array(landmark.loc[3 + (i * 4), ['XM', 'YM']])

                # get rotated bead position
                self._bead_rotated_position[self._patch == patch] = rotate_position(position, target_start, target_end,
                                                                                    initial_start, initial_end)

        else:
            raise ValueError(f'number of rows in the register reference file is nor correct. '
                             f'should be either 4 lines, or 4*n_patch lines.')

    def _link_bead(self):
        print('register bead to well...')

        # register beads to wells, size=well_num, value=bead_index
        d_position, d_inrange = register(self._bead_rotated_position, self._well_position, self._d_th)
        self._well_bead_link = np.repeat(-1, self._well_num)
        self._well_bead_link[d_position[d_inrange]] = np.arange(d_position.size)[d_inrange]

    def _plt_register(self):
        pdf_outfile = self._analysis_folder + 'bead_cell_register.pdf'
        with PdfPages(pdf_outfile) as pdf:
            plt.scatter(self._bead_rotated_position['XM'], self._bead_rotated_position['YM'], c='b', s=0.2, linewidth=0)
            plt.scatter(self._well_position['XM'], self._well_position['YM'], c='r', s=0.4, linewidth=0)
            pdf.savefig()
            plt.close()

    def link_obc_cell(self):
        print("linking bead optical barcode to cell...")

        self._check_params()

        # map bead demultiplexing scan to well scan
        self._rotate_bead_position()

        # register beads to wells
        self._link_bead()

        # check registration
        self._plt_register()

        # get wells with one bead registered, and one cell registered
        values, index, counts = np.unique(self._well_bead_link, return_index=True, return_counts=True)
        well_w_one_bead = set(index[counts == 1])
        well_w_one_cell = set(self.cell.index)
        well_w_bead_cell = set.intersection(well_w_one_bead, well_w_one_cell)
        well_w_bead_cell = list(well_w_bead_cell)

        print('register optical bead barcode to cell...')

        self.obc_cell = pd.DataFrame(-1, columns=['bead_index', 'obc', 'cell_index'],
                                     index=range(len(well_w_bead_cell)))

        # get bead index
        self.obc_cell['bead_index'] = np.array(self._well_bead_link)[well_w_bead_cell]
        # get obc
        self.obc_cell['obc'] = np.array(self._obc.loc[self.obc_cell['bead_index']])
        # get cell_index, well_index==cell_index
        self.obc_cell['cell_index'] = well_w_bead_cell

        # get info
        cell_num = len(well_w_bead_cell)
        print(f'{cell_num} wells have decoded bead-cell pair.')
