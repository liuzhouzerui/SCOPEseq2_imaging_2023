import os
import fnmatch
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from scopeseq.method import register, find_rotation_matrix


class CellIntensity:
    def __init__(self, project_folder):
        self._project_folder = project_folder
        self._bead_folder = None
        self._well_folder = project_folder + 'well/'
        self._cell_folder = None
        self._single_cell_folder = None
        self._analysis_folder = project_folder + 'analysis/'
        self._well_channel = None
        self.n_lane = None
        self.total_lanes = None
        self.d_th = None
        self._well = {}
        self.well_num = None
        self.well_position = None
        self._well_filtered = None
        self._cell_channel = None
        self.cell = None

    def set_project_folder(self, project_folder):
        self._project_folder = project_folder
        self._cell_folder = project_folder+'cell/'
        self._single_cell_folder = project_folder+'cell/single_cell/'
        self._analysis_folder = project_folder + 'analysis/'

    def set_device(self, n_lane=0, total_lanes=1, d_th=40):
        """

        :param n_lane: int - index of device lane
        :param total_lanes: int - number of lanes in the device. 1 for small, medium small device. 5 for large device.
        :param d_th: int - radius of well (pixel)
        :return:
        """
        self.n_lane = n_lane
        self.total_lanes = total_lanes
        self.d_th = d_th

    def set_well_scan(self, well_channel=None):
        """

        :param well_channel: list - fluorescent channels used in cell scan
        :return:
        """
        if well_channel is None:
            self._well_channel = ['GFP', 'TRITC', 'CY5']
        else:
            self._well_channel = well_channel

    def _check_params(self):
        # folders
        print('project_folder: ' + self._project_folder)
        if not os.path.exists(self._project_folder):
            raise ValueError(f'project_folder does not exist')

        print('well_folder: ' + self._well_folder)
        if not os.path.exists(self._well_folder):
            raise ValueError(f'well_folder does not exist')
        print('well_folder should contain '
              '(1) well intensity files from cell scan')

        print('analysis_folder: ' + self._analysis_folder)
        if not os.path.exists(self._analysis_folder):
            os.mkdir(self._analysis_folder)
            print('make analysis folder fot output at ' + self._analysis_folder)

        # values
        if (self.n_lane != int(self.n_lane)) | (self.n_lane < 0) | (
                self.total_lanes != int(self.total_lanes)) | (self.total_lanes < 0):
            raise ValueError(f'n_lane, total_lanes should be non-negative integer.')

        if self.d_th < 0:
            raise ValueError(f'd_th should be non-negative value.')

        if not isinstance(self._well_channel, list):
            raise TypeError(f'well_channel should be a list.')

    def _find_well_intensity_file(self, channel):
        # generate file pattern
        file_seq = str(10000 + self.n_lane)[1:5]
        file_ext = channel + ".tif_Results.xls"
        file_pattern = '*' + file_seq + '*' + file_ext

        # find file
        file_name_list = fnmatch.filter(os.listdir(self._well_folder), file_pattern)
        if len(file_name_list) == 0:
            raise OSError(f'The well intensity file is not found!')
        else:
            file_name = file_name_list[0]

        return file_name

    def _initialize(self, channel):
        """

        :param channel: channel name used for the well image
        :return:
        """
        print("initializing...")

        # get initialize intensity matrix
        file_name = self._find_well_intensity_file(channel)
        print(file_name)
        intensity_matrix = pd.read_csv(self._well_folder + file_name, sep="\t")

        # get number of wells
        self.well_num = intensity_matrix.shape[0]

        # get well position
        self.well_position = intensity_matrix[['XM', 'YM']]

        # delete intensity matrix
        del intensity_matrix

    def _register(self, channel):
        print("generating well-based cell intensity..." + channel)

        # get intensity matrix
        file_name = self._find_well_intensity_file(channel)
        print(file_name)
        intensity_matrix = pd.read_csv(self._well_folder + file_name, sep="\t")

        # record
        self._well[channel] = intensity_matrix

        # delete intensity matrix
        del intensity_matrix

    def _plt_intensity(self):
        pdf_outfile = self._analysis_folder + 'well_intensity.pdf'
        with PdfPages(pdf_outfile) as pdf:
            for c in self._well_channel:
                plt.hist(np.log2(self._well[c]['Mean']), bins=200)
                plt.title(c)
                pdf.savefig()
                plt.close()

    def generate_well_intensity(self):
        self._check_params()

        # intialization
        self._initialize(self._well_channel[0])

        # generate well intensity
        for c in self._well_channel:
            self._register(c)

        # plt well intensity
        self._plt_intensity()

    def cell_filter(self, filter_th):
        """

        :param filter_th: dict - {fluorescent channel: threshold}
        :return:
        """
        # check parameter
        if not isinstance(filter_th, dict):
            raise TypeError(f'filter_th should be a dictionary.')

        for c in filter_th.keys():
            if c not in self._well_channel:
                raise ValueError(f'keys in filter_th should match values in well_channel.')

        # get well intensity
        well_header = []
        for c in self._well.keys():
            well_header.append(c)
        well = pd.DataFrame(None, index=range(self.well_num), columns=well_header)

        for c in self._well.keys():
            well[c] = self._well[c]['Mean']

        # get well information
        well['index'] = range(self.well_num)
        well['well_x'] = self.well_position['XM']
        well['well_y'] = self.well_position['YM']

        # filter cell
        well_have_cells = set()
        for c in filter_th.keys():
            well_have_cells = well_have_cells.union(well['index'][well[c].values > filter_th[c]])
        self._well_filtered = well.loc[sorted(list(well_have_cells))]

        # output
        self._well_filtered.to_csv(self._well_folder + 'Results.csv', index=False)

        # get info
        cell_num = len(well_have_cells)
        print(f'{cell_num} wells have cells.')

    def generate_cell_intensity(self, cell_channel):
        """

        :param cell_channel: list - fluorescent channels used in cell scan
        :return:
        """
        # get cell measurement path
        self._cell_folder = self._project_folder + 'cell/'
        self._single_cell_folder = self._cell_folder + 'single_cell/'

        self._cell_channel = cell_channel

        # folders
        print('cell_folder: ' + self._cell_folder)
        if not os.path.exists(self._cell_folder):
            raise ValueError(f'cell_folder does not exist')

        for c in self._cell_channel:
            path = self._cell_folder + c + '/intensity_data/'
            print('cell_folder_' + c + ': ' + path)
            if not os.path.exists(path):
                raise ValueError(f'cell_folder_{c} does not exist')

        print('single_cell_folder: ' + self._single_cell_folder)
        if not os.path.exists(self._single_cell_folder):
            raise ValueError(f'single_cell_folder does not exist')

        # get measurement parameters
        index = self._well_filtered.index
        columns = self._well[list(self._well.keys())[0]].columns

        # get cell_based single cell intensity
        image_features = {}
        for c in self._cell_channel:
            print('generate cell intensity...' + c)

            # initialize
            image_features[c] = pd.DataFrame(None, index=index, columns=columns)
            image_features[c + '_bg'] = pd.DataFrame(None, index=index, columns=columns)

            # get measurement
            for i in index:
                file = self._cell_folder + c + '/intensity_data/' + str(i) + '.' + c + '.tif_Results.xls'
                if os.path.exists(file):
                    if os.path.getsize(file) > 0:
                        intensity_matrix = pd.read_csv(file, sep='\t')
                        if intensity_matrix.shape[0] == 2:
                            image_features[c].loc[i, columns] = intensity_matrix.loc[0, columns]
                            image_features[c + '_bg'].loc[i, columns] = intensity_matrix.loc[1, columns]

        # rearrange cell intensity
        self.cell = pd.DataFrame(None, index=index)
        for c in self._cell_channel:
            for column in columns:
                self.cell[c + '_' + column] = image_features[c][column].values
        for c in self._cell_channel:
            for column in columns:
                self.cell[c + '_bg_' + column] = image_features[c + '_bg'][column].values
            self.cell[c + '_Mean_norm'] = self.cell[c + '_Mean'] - self.cell[c + '_bg_Mean']

        # single cell filter
        # image filter
        single_cell = os.listdir(self._single_cell_folder)
        if '.DS_Store' in single_cell:
            single_cell.remove('.DS_Store')
        single_cell = [x.split('.')[0] for x in single_cell]
        self.cell = self.cell[np.isin(index, single_cell)]
        # ImageJ filter
        for c in self._cell_channel:
            self.cell = self.cell[self.cell[c + '_Mean_norm'].notna()]

        self.cell.to_csv(self._cell_folder + 'sc_image.matrix.txt', sep='\t')

        cell_num = self.cell.shape[0]
        print(f'{cell_num} single cells.')
