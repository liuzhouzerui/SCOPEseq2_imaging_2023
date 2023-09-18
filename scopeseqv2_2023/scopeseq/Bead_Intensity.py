import os
import fnmatch
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from scopeseq.method import register, assign_obc, assign_obc_stage, quantile_linear_transformation, obc_correlation
from scopeseq.stats import cutoff, fit_double_gaussian


class BeadIntensity:
    def __init__(self, project_folder):
        self._project_folder = project_folder
        self._bead_folder = project_folder+'bead/intensity_data/'
        self._analysis_folder = project_folder + 'analysis/'
        self.n_lane = None
        self.total_lanes = None
        self.n_patch_per_lane = None
        self._d_th = None
        self._first_border = None
        self._border_gap = None
        self._bead_channel = None
        self._round = None
        self._barcode_ref_fn = None
        self.bead_num = None
        self._bead_area = None
        self.bead_position = None
        self._probe = None
        self._probe_raw = None
        self._probe_normalized = None
        self._bg = None
        self.patch = np.array([])
        self._reverse = False
        self._obc_calling = None
        self.obc = []
        self._obc_s = []
        self._obc_q = []
        self._obc_round = None
        self._correlation_s = None
        self._correlation_q = None
        self._obc_coef_s = None
        self._obc_coef_q = None

    def set_project_folder(self, project_folder, bead_folder=None, analysis_folder=None):
        self._project_folder = project_folder
        if bead_folder is not None:
            self._bead_folder = bead_folder
        if bead_folder is not None:
            self._analysis_folder = analysis_folder

    def set_device(self, n_lane=0, total_lanes=1, patches_per_lane=10, first_border=600, border_gap=7400, d_th=40,
                   reverse=False):
        """

        :param n_lane: int - index of device lane
        :param total_lanes: int - number of lanes in the device. 1 for small, medium small device. 5 for large device.
        :param patches_per_lane: int - number of segments per device lane. 10 for medium small. 6 for large segment.
        :param first_border: int - x pixel position for the start of the first segment.
        :param border_gap: int - segment length in pixel along the x axis
        :param d_th: int - radius of well (pixel) - radius of bead (pixel)
        :return:
        """
        self.n_lane = n_lane
        self.total_lanes = total_lanes
        self.n_patch_per_lane = patches_per_lane
        self._first_border = first_border
        self._border_gap = border_gap
        self._d_th = d_th
        self._reverse = reverse

    def set_demultiplexing(self, round_number=8, bead_channel=None, barcode_ref_fn=None):
        """

        :param round_number: int - total number of demultiplexing cycles
        :param bead_channel: list - fluorescent channels used in demultiplexing
        :param barcode_ref_fn: file name - bead optical barcode, sequence to binary code reference
        :return:
        """
        # get bead imaging channel names
        if bead_channel is None:
            self._bead_channel = ['bf', 'cy5', 'cy3']
        else:
            self._bead_channel = bead_channel

        # get demultiplexing cycles
        self._round = round_number

        # get optical barcode reference
        if barcode_ref_fn is None:
            self._barcode_ref_fn = self._bead_folder + '192_8mer_seq_reference.csv'
        else:
            self._barcode_ref_fn = barcode_ref_fn

    def _check_params(self):
        # folders
        print('project_folder: ' + self._project_folder)
        if not os.path.exists(self._project_folder):
            raise ValueError(f'project_folder does not exist')

        print('bead_folder: ' + self._bead_folder)
        if not os.path.exists(self._bead_folder):
            raise ValueError(f'bead_folder does not exist')
        print('bead_folder should contain '
              '(1) bead intensity files from optical demultiplexing scans '
              '(2) barcode reference')

        print('analysis_folder: ' + self._analysis_folder)
        if not os.path.exists(self._analysis_folder):
            os.mkdir(self._analysis_folder)
            print('make analysis folder fot output at ' + self._analysis_folder)

        # file
        print('barcode_ref_fn: ' + self._barcode_ref_fn)
        if not os.path.exists(self._barcode_ref_fn):
            raise ValueError(f'barcode_ref_fn does not exist')

        # values
        if (self.n_lane != int(self.n_lane)) | (self.n_lane < 0) | (
                self.total_lanes != int(self.total_lanes)) | (self.total_lanes < 0) | (
                self.n_patch_per_lane != int(self.n_patch_per_lane)) | (self.n_patch_per_lane < 0):
            raise ValueError(f'n_lane, total_lanes, patches_per_lane should be non-negative integer.')

        if (self._first_border < 0) | (
                self._border_gap < 0) | (
                self._d_th < 0):
            raise ValueError(f'first_border, border_gap, d_th should be non-negative value.')

        if (not isinstance(self._bead_channel, list)) | (len(self._bead_channel) != 3):
            raise ValueError(f'bead_channel should be a three value list. '
                             f'ordered as bright field, cy5 (S probe), cy3 (Q probe)')

        if (self._round != int(self._round)) | (self._round < 0):
            raise ValueError(f'round should be non-negative integer.')

    def _find_bead_intensity_file(self, iter_round, channel, bp):
        """

        :param iter_round:
        :param channel:
        :param bp:
        :return:
        """
        # generate file pattern
        file_round = '0000' + str(iter_round - 1)
        if bp == 'bg':
            file_seq = str(10000 + (iter_round - 1) * self.total_lanes * 2 + self.n_lane)[1:5]
        elif bp == 'probe':
            file_seq = str(10000 + (iter_round - 1) * self.total_lanes * 2 + self.total_lanes + self.n_lane)[1:5]
        else:
            file_seq = ''
        file_ext = self._bead_channel[channel] + ".tif_Results.xls"
        file_pattern = '*' + file_round + '*' + file_seq + '*' + file_ext

        # find file
        file_name_list = fnmatch.filter(os.listdir(self._bead_folder), file_pattern)
        if len(file_name_list) == 0:
            raise OSError(f'The bead intensity file is not found!')
        else:
            file_name = file_name_list[0]

        return file_name

    def _initialize(self, iter_round, channel, bp):
        """

        :param iter_round:
        :param channel:
        :param bp:
        :return:
        """
        print("initializing...")

        # get initialize intensity matrix
        file_name = self._find_bead_intensity_file(iter_round, channel, bp)
        print(file_name)
        intensity_matrix = pd.read_csv(self._bead_folder + file_name, sep='\t')

        # get bead number
        self.bead_num = intensity_matrix.shape[0]

        # get bead area
        self._bead_area = intensity_matrix['Area']

        # get bead position
        self.bead_position = intensity_matrix[['XM', 'YM']]

        # initialize probe intensity matrix
        probe_header = []
        for c in self._bead_channel[1:3]:
            for i in range(1, self._round+1, 1):
                probe_header.append('probe_' + str(i) + '_' + c)
        self._probe = pd.DataFrame(-1, columns=probe_header, index=range(self.bead_num))

        # initialize background intensity matrix
        bg_header = []
        for c in self._bead_channel[1:3]:
            for i in range(1, self._round+1, 1):
                bg_header.append('bg_' + str(i) + '_' + c)
        self._bg = pd.DataFrame(-1, columns=bg_header, index=range(self.bead_num))

        # assign intensity
        if bp == 'bg':
            self._bg[bp + '_' + str(iter_round) + '_' + self._bead_channel[channel]] = \
                intensity_matrix['Mean'].values
        if bp == 'probe':
            self._probe[bp + '_' + str(iter_round) + '_' + self._bead_channel[channel]] = \
                intensity_matrix['Mean'].values

        # delete intensity matrix
        del intensity_matrix

    def _register_bead(self, iter_round, channel, bp):
        """

        :param iter_round:
        :param channel:
        :param bp:
        :return:
        """
        print("registering..." + bp + '_' + str(iter_round) + '_' + self._bead_channel[channel])

        # get following intensity matrix
        file_name = self._find_bead_intensity_file(iter_round, channel, bp)
        print(file_name)
        intensity_matrix = pd.read_csv(self._bead_folder + file_name, sep='\t')

        # register to the first bead demultiplexing image
        d_position, d_inrange = register(intensity_matrix[['XM', 'YM']], self.bead_position, self._d_th)

        # assign intensity
        if bp == 'bg':
            self._bg.loc[d_position[d_inrange], bp + '_' + str(iter_round) + '_' + self._bead_channel[channel]] \
                = intensity_matrix.iloc[np.arange(d_position.size)[d_inrange], 2].values
        if bp == 'probe':
            self._probe.loc[d_position[d_inrange], bp + '_' + str(iter_round) + '_' + self._bead_channel[channel]] \
                = intensity_matrix.iloc[np.arange(d_position.size)[d_inrange], 2].values

        # delete intensity matrix
        del intensity_matrix

    def _assign_patch(self):
        print('assigning patches...')
        borders = np.arange(self._first_border, self._first_border + (self.n_patch_per_lane + 1) * self._border_gap,
                            self._border_gap)
        patch = self.bead_position['XM'].apply(lambda x: np.argmax((borders > x) * 1)-1)
        if not self._reverse:
            self.patch = patch + self.n_patch_per_lane * self.n_lane
        if self._reverse:
            self.patch = self.n_patch_per_lane * (self.n_lane+1) - patch -1

        # plt patch information to check
        pdf_outfile = self._analysis_folder+'bead_patch.pdf'
        with PdfPages(pdf_outfile) as pdf:
            plt.scatter(self.bead_position['XM'], self.bead_position['YM'], c=self.patch, s=0.5)
            for x in borders:
                plt.vlines(x, ymin=self.bead_position['YM'].min(), ymax=self.bead_position['YM'].max())
            pdf.savefig()
            plt.close()

    def _plt_intensity(self):
        """
        plot bead intensities by cycle
        :return:
        """
        pdf_outfile = self._analysis_folder + 'bead_intensity.pdf'
        with PdfPages(pdf_outfile) as pdf:
            for i in range(self._probe.shape[1]):
                plt.hist(np.log2(self._probe.iloc[:, i][self._probe.iloc[:, i] > 0]), bins=50)
                plt.title(self._probe.columns[i])
                pdf.savefig()
                plt.close()

    def generate_bead_intensity(self):
        """
        de-multiplexing image intensity.
        :return:
        """
        self._check_params()

        for iter_round in range(1, self._round+1, 1):
            for channel in [1, 2]:
                for bp in ['bg', 'probe']:
                    if self.bead_position is None:
                        self._initialize(iter_round, channel, bp)
                    else:
                        self._register_bead(iter_round, channel, bp)

        # assign patch
        self._assign_patch()

        # get plots
        self._plt_intensity()

        # get info
        bead_num = sum(self._probe.min(axis=1) > 0)
        print(f'get {self.bead_num} beads, {bead_num} successfully registered.')

    def _intensity_normalization(self, method='linear', no_signal_th=1000):
        if method not in ['ql', 'linear']:
            raise ValueError(f'bead intensity normalization method should be "linear".')

        if method == 'linear':
            probe = self._probe
            probe[probe.max(axis=1) < no_signal_th] = -1
            self._probe_normalized = probe.apply(lambda x: x if x.min() == -1 else (x - x.min()) / x.max(), axis=1)

    def probe_norm(self):
        """
        normalize probe intensities across cycles
        :return:
        """
        self._probe_raw = self._probe.copy()

        bg = {}
        for i in self._probe.columns[0:self._round]:
            x = self._probe_raw[i]
            x = np.log2(x[(~np.isnan(x)) & (x > 1)])
            chist, cedge = np.histogram(x, bins=50)
            md = np.median(x)
            mn1 = cedge[np.argmax(chist[cedge[0:-1] < md])]
            bg[i] = mn1
        norm = min(bg.values())
        for i in self._probe.columns[0:self._round]:
            self._probe[i] = self._probe_raw[i].apply(lambda x: x/(2**(bg[i]-norm)))

        bg = {}
        for i in self._probe.columns[self._round:self._probe.shape[1]]:
            x = self._probe_raw[i]
            x = np.log2(x[(~np.isnan(x)) & (x > 1)])
            chist, cedge = np.histogram(x, bins=50)
            md = np.median(x)
            mn1 = cedge[np.argmax(chist[cedge[0:-1] < md])]
            bg[i] = mn1
        norm = min(bg.values())
        for i in self._probe.columns[self._round:self._probe.shape[1]]:
            self._probe[i] = self._probe_raw[i].apply(lambda x: x/(2**(bg[i]-norm)))

        self._plt_intensity()

    def obc_calling_core(self, no_signal_th=None, mode='all', normalize=None):
        """
        barcode calling using core-by-core method (bead-by-bead)
        :param no_signal_th: int - the minimum intensity required for demultiplexing
        :param mode: category - 'max' is one round core-by-core method, 'all' is full round core-by-core method
        :param normalize: category -  'ql' or 'linear'
        :return:
        """
        print('optical barcode calling using core-by-core method...')
        self._obc_calling = 'core'

        # get barcode reference
        barcode_ref_table = pd.read_csv(self._barcode_ref_fn, dtype=str)

        # probe intensity normalization
        if normalize is None:
            probe = self._probe
        else:
            self._intensity_normalization(method=normalize)
            probe = self._probe_normalized

        # demultiplex s probe
        self._obc_s = probe.iloc[:, 0:self._round].apply(
            lambda x: assign_obc(x, barcode_ref_table['Barcode_S'], no_signal_th=no_signal_th, mode=mode),
            axis=1, result_type="expand")

        # demultiplex q probe
        self._obc_q = probe.iloc[:, self._round:self._probe.shape[1]].apply(
            lambda x: assign_obc(x, barcode_ref_table[
                'Barcode_Q'], no_signal_th=no_signal_th, mode=mode), axis=1, result_type="expand")

        # generate bead optical barcode
        self.obc = self.patch.astype('str') + '_' + self._obc_s[0].astype('str') + '_' + self._obc_q[0].astype('str')
        self._obc_round = self._obc_s[1].astype('str') + '_' + self._obc_q[1].astype('str')

        # get plots
        self._plt_obc()

    def obc_calling_stage(self):
        """
        barcode calling using stage-by-stage method (cycle-by-cycle)
        :return:
        """
        print('optical barcode calling using stage-by-stage method...')
        self._obc_calling = 'stage'

        # get barcode reference
        barcode_ref_table = pd.read_csv(self._barcode_ref_fn, dtype=str)

        obc = pd.DataFrame(None, index=self._probe.index, columns=self._probe.columns)

        # determine cut-off of each stage and assign on-off
        for i in self._probe.columns:
            probe = self._probe[i]
            log_probe = np.log2(probe)
            # get cutoff threshold
            thresh = cutoff(probe)
            # assign on-off state
            obc[i][(log_probe > 0) & (log_probe <= thresh)] = 0
            obc[i][log_probe > thresh] = 1
            obc[i][probe == -1] = -1

        # demultiplex s probe
        self._obc_s = obc.iloc[:, 0:self._round].apply(
            lambda x: assign_obc_stage(x, barcode_ref_table['Barcode_S']), axis=1).values

        # demultiplex q probe
        self._obc_q = obc.iloc[:, self._round:self._probe.shape[1]].apply(
            lambda x: assign_obc_stage(x, barcode_ref_table['Barcode_Q']), axis=1).values

        # get bead optical barcode
        self.obc = self.patch.astype('str') + '_' + self._obc_s_stage.astype(
            'str') + '_' + self._obc_q_stage.astype('str')
        self._obc_round = self._obc

        # get plots
        self._plt_obc()

    def obc_calling_correlation(self, no_signal_th=1000):
        """
        barcode calling using pearson correlation
        :return:
        """
        print('optical barcode calling using correlation method...')
        self._obc_calling = 'correlation'

        # get barcode reference
        barcode_ref_table = pd.read_csv(self._barcode_ref_fn, dtype=str)

        ref_header = []
        for i in range(1, self._round + 1, 1):
            ref_header.append('Barcode_S_' + str(i))
        barcode_ref_s = barcode_ref_table[ref_header].astype('int')
        ref_header = []
        for i in range(1, self._round + 1, 1):
            ref_header.append('Barcode_Q_' + str(i))
        barcode_ref_q = barcode_ref_table[ref_header].astype('int')

        # probe intensity normalization
        print('normalize intensity...')
        self._intensity_normalization(method='linear', no_signal_th=no_signal_th)

        # calculate correlation for barcode s
        print('calculate correlation...')
        self._correlation_s = self._probe_normalized.iloc[:, 0:self._round].apply(
            lambda x: obc_correlation(x, barcode_ref_s), axis=1)
        # assign barcode by the highest Pearson coefficient
        self._obc_s = self._correlation_s.apply(lambda x: -1 if x.max() == -1 else x.values.argmax(), axis=1)
        # record the the highest Pearson coefficient as assignment confidence
        self._obc_coef_s = self._correlation_s.apply(lambda x: -1 if x.max() == -1 else x.values.max(), axis=1)

        # calculate correlation for barcode q
        self._correlation_q = self._probe_normalized.iloc[:, self._round:self._probe_normalized.shape[1]].apply(
            lambda x: obc_correlation(x, barcode_ref_q), axis=1)
        self._obc_q = self._correlation_q.apply(lambda x: -1 if x.max() == -1 else x.values.argmax(), axis=1)
        self._obc_coef_q = self._correlation_q.apply(lambda x: -1 if x.max() == -1 else x.values.max(), axis=1)

        # generate bead optical barcode
        self.obc = self.patch.astype('str') + '_' + self._obc_s.astype('str') + '_' + self._obc_q.astype('str')
        self._obc_round = self.obc

        # get plots
        self._plt_obc()

    def _plt_obc(self):
        if self._obc_s.size == self.bead_num:
            s = self._obc_s
            q = self._obc_q
        else:
            s = self._obc_s[0]
            q = self._obc_q[0]

        pdf_outfile = self._analysis_folder + 'obc_usage.pdf'
        with PdfPages(pdf_outfile) as pdf:
            plt.hist(s, bins=np.arange(0, 97, 1))
            plt.title('S Barcode')
            pdf.savefig()
            plt.close()

            plt.hist(q, bins=np.arange(0, 97, 1))
            plt.title('Q Barcode')
            pdf.savefig()
            plt.close()
