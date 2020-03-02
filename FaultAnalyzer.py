import pandas as pd
import numpy as np
import scipy as sp
import scipy.signal as sg
from scipy.fftpack import fft as spfft
from scipy.fftpack import fftfreq as spfftfreq
import matplotlib as mp
import matplotlib.pyplot as plt

# Following code only needed for latex compatible plots
# size = 20
# pgf_with_latex = {                      # setup matplotlib to use latex for output
#     "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
#     "text.usetex": True,                # use LaTeX to write all text
#     "font.family": 'serif',
#     "font.serif": [],                   # blank entries should cause plots
#     "font.sans-serif": [],              # to inherit fonts from the document
#     "font.monospace": [],
#     "axes.labelsize": size,               # LaTeX default is 10pt font.
#     "font.size": size,
#     "legend.fontsize": size,               # Make the legend/label fonts
#     "xtick.labelsize": size,               # a little smaller
#     "ytick.labelsize": size,
#     # "figure.figsize": (12, 8),     # default fig size of 0.9 textwidth
#     "pgf.preamble": [
#         r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts
#         r"\usepackage[T1]{fontenc}",        # plots will be generated
#         r"\usepackage[detect-all,locale=DE]{siunitx}",
#         ]                                   # using this preamble
#     }
# # }}}
# mp.rcParams.update(pgf_with_latex)


def convert_data():
    """Convert raw .mat data into compressed pickle."""
    from raw_data.mat_to_pickle import to_pickle
    to_pickle()


class Analyzer:
    def __init__(self, run_park_tr=False, run_fft=False, fault_display=['stf'], upper_freq_lim=500):
        self.data = {
            'h': pd.read_pickle('data_healthy.pickle'),
            'f1': pd.read_pickle('fault1.pickle'),
            'f2': pd.read_pickle('fault2.pickle'),
            'f3': pd.read_pickle('fault3.pickle'),
        }
        self.data['f1']['t'] -= self.data['f1']['t'][0]
        self.numpy_data = {}
        for key in self.data.keys():
            self.numpy_data[key] = self.data[key].transpose().iloc[1:].to_numpy()

        # initial assignments
        self.def_data_names = ['h', 'f1', 'f2', 'f3']
        self.freq_h1 = {
            key: 48.8 for key in self.def_data_names
        }
        self.magn_h1 = {}
        self.freq_mech = 1420/60
        self.def_data_names_col = {
            'h': 'b',
            'f1': 'orange',
            'f2': 'g',
            'f3': 'r',
        }
        self.upper_freq_lim = upper_freq_lim
        self.time = {}      # time DataFrames
        for key in self.def_data_names:
            self.time[key] = self.data[key]['t']
        self.data_len = {}  # data lengths
        for key in self.def_data_names:
            self.data_len[key] = self.time[key].count()
        self.data_sp = {}   # data sample times
        for key in self.def_data_names:
            self.data_sp[key] = self.time[key][1]
        self.park_data = {}
        self.numpy_park_data = {}
        self.fft_data = {}
        self._process_input(run_park_tr, run_fft, fault_display)

    def _process_input(self, run_park_tr, run_fft, fault_display):
        """Process data if needed."""
        if run_park_tr:
            park_data = self._run_extended_park_transformation()
            for key in self.def_data_names:
                self.park_data[key] = pd.DataFrame(park_data[key])
                self.park_data[key].to_pickle(str('park_' + key + '.pickle'))
                self.numpy_park_data[key] = self.park_data[key].transpose().to_numpy()
        elif run_fft:
            for key in self.def_data_names:
                self.park_data[key] = pd.read_pickle(str('park_' + key + '.pickle'))
                self.numpy_park_data[key] = self.park_data[key].transpose().to_numpy()
        if run_fft:
            fft_data = self._run_fft()
            for key in self.def_data_names:
                self.fft_data[key] = pd.DataFrame(fft_data[key])
                self.fft_data[key].to_pickle(str('fft_' + key + '.pickle'))
        else:
            for key in self.def_data_names:
                self.fft_data[key] = pd.read_pickle(str('fft_' + key + '.pickle'))
        for key in self.def_data_names:
            self.fft_data[key] = self.fft_data[key].loc[
                (self.fft_data[key]['freq'] <= self.upper_freq_lim)
            ].reset_index(drop=True)
        self._find_fundamental_el_freq(self.freq_h1['h'])
        if fault_display:
            self.fault_freq, self.fault_freq_style = self._fault_freq_detection(fault_display)
        self.magn_faults = {}
        self.freq_faults = {}
        self._find_significant_fault_val(fault_display)

    def _run_extended_park_transformation(self):
        """Run extended park transformation."""
        data = self.numpy_data
        park_data = {}
        for key in self.def_data_names:
            i_d = np.sqrt(2 / 3) * data[key][0] - np.sqrt(1 / 6) * data[key][1] - np.sqrt(1 / 6) * data[key][2]
            i_q = np.sqrt(1 / 2) * data[key][1] - np.sqrt(1 / 2) * data[key][2]
            i_p = np.sqrt(i_d ** 2 + i_q ** 2)
            park_data[key] = {
                'id': i_d,
                'iq': i_q,
                'ip': i_p,
            }
        return park_data

    def _run_fft(self):
        """Run fft algorithm on data

        Optionally run when creating FaultAnalyzer.
        Additionally calculating base frequency for each data series.
        Rewriting 'fft_data.pickle'.
        """
        data = self.numpy_park_data
        fft_data = {}
        for key in self.def_data_names:
            tmp = data[key]
            le = self.data_len[key]
            if le % 2:
                le -= 1
            fft_data[key] = {}
            for idx in range(3):
                fft = 2 * spfft(tmp[idx], le) / le
                fft = abs(fft)
                fft = fft[range(int(le / 2))]
                key2 = list(self.park_data[key].keys())
                fft_data[key][key2[idx]] = fft
            fft_data[key]['freq'] = spfftfreq(self.data_len[key], d=self.data_sp[key])[range(int(le / 2))]
        return fft_data

    def _find_fundamental_el_freq(self, start_freq):
        """Computing fundamental frequencies from data."""
        range_freq = 5
        for key in self.fft_data:
            mx, freq = self._find_peak_freq(
                start_freq,
                range_freq,
                self.fft_data[key],
            )
            self.freq_h1[key] = np.mean(list(freq.values()))
            self.magn_h1[key] = mx

    def _find_peak_freq(self, start_freq, range_freq, data):
        """
        Calculate base frequency of data within given frequency range of rated electric frequency
        """
        tmp_data = data.loc[
            data['freq'].between(start_freq - range_freq, start_freq + range_freq)
        ]
        tmp = tmp_data.filter(regex='i')[2:].idxmax()
        mx = {
            'id': tmp_data['id'][tmp[0]],
            'iq': tmp_data['iq'][tmp[1]],
            'ip': tmp_data['ip'][tmp[2]],
        }
        tmp = [tmp_data['freq'].loc[item] for item in tmp]
        freq = {
            'id': tmp[0],
            'iq': tmp[1],
            'ip': tmp[2],
        }
        return mx, freq

    def _find_significant_fault_val(self, fault_display):
        """Compute magnitudes around calculated fault frequencies.

        Parameters
        ----------
        fault_display: list of str
               Names of faults to display
        """
        range_freq = 3
        magn_fault = {}
        freq_fault = {}
        for f in fault_display:
            faults = self.fault_freq[f]
            magn_fault[f] = {}
            freq_fault[f] = {}
            for item in ['id', 'iq', 'ip']:
                magn_fault[f][item] = {}
                freq_fault[f][item] = {}
                for key in self.def_data_names:
                    magn_fault[f][item][key] = []
                    freq_fault[f][item][key] = []
                    for tmp in faults[key][:-1]:
                        mx, freq = self._find_peak_freq(tmp, range_freq, self.fft_data[key])
                        magn_fault[f][item][key].append(mx[item])
                        freq_fault[f][item][key].append(freq[item])
        self.magn_faults = magn_fault
        self.freq_faults = freq_fault

    def _fault_freq_detection(self, fault_display):
        """Evaluate fault frequencies.

        Parameters
        ----------
        fault_display: list of str
               Names of faults to display
        """
        faults = {
            'bearing': [],  # self._fault_bearing,
            'stf': self._fault_stf,
        }
        fault_freq = {}
        fault_freq_style = {}
        tmp_style = ['-.', ':', '-.-', '--']
        idx = 0
        for fault in fault_display:
            fault_freq[fault] = faults[fault]()
            fault_freq_style[fault] = tmp_style[idx]
            idx += 1
        return fault_freq, fault_freq_style

    def _fault_stf(self):
        """Short turn fault freq based on Chapter4/Slide20"""
        tmp = {}
        for key in self.def_data_names:
            i = 0
            val = 0
            tmp[key] = []
            while val < self.upper_freq_lim:
                val = self.freq_h1[key] * i
                tmp[key].append(val)
                i += 1
        return tmp

    def plot_stft(self, item=2):
        """Plot stft of park components of all data sets."""
        # loc = '../sigproc_report/'
        fig = plt.figure(figsize=(24, 12))
        gs = mp.gridspec.GridSpec(2, 3, figure=fig)
        gs.set_width_ratios([2, 2, .2])
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[0, 1])
        ax4 = fig.add_subplot(gs[1, 1])
        ax5 = fig.add_subplot(gs[:, 2])
        axes = {'h': ax1, 'f1': ax2, 'f2': ax3, 'f3': ax4}
        data = self.numpy_park_data
        fault_names = {
            'h': 'H',
            'f1': 'F1',
            'f2': 'F2',
            'f3': 'F3',
        }
        components = {
            0: '$i_d$',
            1: '$i_q$',
            2: '$i_p$',
        }
        for key in data.keys():
            ax = axes[key]
            fs = self.data_sp[key] ** -1
            f, t, zxx = sg.stft(data[key][item], fs, nperseg=40000)
            zxx = np.abs(zxx)
            im = ax.pcolormesh(
                t,
                f,
                zxx,
                vmin=.0001,
                vmax=.007,
                cmap='copper_r',
                edgecolors='None',
            )

            ax.title.set_text(fault_names[key])
            ax.set_ylim(0, 200)
            ax.set_yticks([0, 50, 100, 150, 200])
            if key == 'f3':
                ax.set_xlabel('time (s)')
            if key == 'h':
                ax.set_ylabel('frequency (hz)')
            if key == 'f1':
                ax.set_ylabel('frequency (hz)')
                ax.set_xlabel('time (s)')
        fig.colorbar(im, cax=ax5)
        ax5.set_ylabel(components[item] + ' amplitude (A)')
        fig.tight_layout()
        plt.show()

    def plot_peak_comparison(self):
        """Plotting magnitudes of fault frequencies."""
        # loc = '../sigproc_report/'
        for f in self.freq_faults.keys():
            fig = plt.figure(figsize=(36, 24))
            gs = mp.gridspec.GridSpec(3, 1, figure=fig)
            ax1 = fig.add_subplot(gs[0, :])
            ax2 = fig.add_subplot(gs[1, :])
            ax3 = fig.add_subplot(gs[2, :])
            components = {
                'id': '$i_d$',
                'iq': '$i_q$',
                'ip': '$i_p$',
            }
            magn = self.magn_faults[f]
            data = {}
            axes = {
                'id': ax1,
                'iq': ax2,
                'ip': ax3,
            }
            for item in magn.keys():
                data[item] = pd.DataFrame(
                    magn[item],
                    index=[str('$h_{' + str(x) + '}$') for x in range(0, len(magn[item]['h']))],
                )
            for item in data.keys():
                ax = axes[item]
                plot_data = data[item].copy()
                plot_data.plot(
                    ax=ax,
                    kind='bar',
                )
                ax.set_ylabel(components[item] + ' amplitude (A)')
                ax.set_xlabel('harmonics')
                names = {
                    'h': 'H',
                    'f1': 'F1',
                    'f2': 'F2',
                    'f3': 'F3',
                }
                ax.legend([val for val in names.values()])
                ax.grid('on')
            fig.tight_layout()
            plt.show()


if __name__ == "__main__":
    """if not existing, create .pickle data from .mat files
    The .mat files have to inserted into 'raw_data' before execution."""
    convert_data()

    """run Analyzer, possible execution of FFT and Park transformation necessary"""
    obj = Analyzer(
        run_park_tr=True,
        run_fft=True,
    )
    """plot STFT, data only displayed, if object created with key 'run_fft=True'
    set 'item=<int>' for component according to: 'ip'=2, 'iq'=1, 'id'=0"""
    # obj.plot_stft(item=1)
    [obj.plot_stft(item=i) for i in [0, 1, 2]]
    """plot peak comparison of extended park vectors"""
    obj.plot_peak_comparison()
