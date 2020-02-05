"""
Created on Mon Jan 22 12:14:10 2020

@author: Ernst
"""

import pandas as pd
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import pickle
import time


class FaultAnalyzer:
    def __init__(
            self,
            run_park_tr=False,
            run_fft=False,
            fault_display=None,
            upper_freq_lim=500
    ):
        self.data = {
            'h': pd.read_pickle('data_healthy.pickle'),
            'f1': pd.read_pickle('fault1.pickle'),
            'f2': pd.read_pickle('fault2.pickle'),
            'f3': pd.read_pickle('fault3.pickle'),
        }
        self.numpy_data = {}
        for key in self.data.keys():
            self.numpy_data[key] = self.data[key].transpose().iloc[1:].to_numpy()

        # initial assignments
        self.def_data_names = ['h', 'f1', 'f2', 'f3']
        self.freq_el = {
            key: 48.8 for key in self.def_data_names
        }
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
        # print('         sample time  | sample length \n' +
        #       'healthy |   ' + str(h_sp) + '     |  ' + str(h_n) + '\n' +
        #       'fault1  |   ' + str(f1_sp) + '     |  ' + str(f1_n) + '\n' +
        #       'fault2  |   ' + str(f2_sp) + '     |  ' + str(f2_n) + '\n' +
        #       'fault3  |   ' + str(f3_sp) + '     |  ' + str(f3_n))

    def _process_input(self, run_park_tr, run_fft, fault_display):
        """Process data if needed."""
        if run_park_tr:
            park_data = self._extended_park_transformation()
            for key in self.def_data_names:
                self.park_data[key] = pd.DataFrame(park_data[key])
                self.park_data[key].to_pickle(str('park_' + key + '.pickle'))
                self.numpy_park_data[key] = self.park_data[key].transpose().to_numpy()
        elif run_fft:
            # with open('park_data.pickle', 'rb') as handle:
            #     self.park_data = pickle.load(handle)
            for key in self.def_data_names:
                self.park_data[key] = pd.read_pickle(str('park_' + key + '.pickle'))
                self.numpy_park_data[key] = self.park_data[key].transpose().to_numpy()
        if run_fft:
            fft_data = self.run_fft()
            for key in self.def_data_names:
                self.fft_data[key] = pd.DataFrame(fft_data[key])
                self.fft_data[key].to_pickle(str('fft_' + key + '.pickle'))
        else:
            # with open('fft_data.pickle', 'rb') as handle:
            #     self.fft_data, self.freq_el = pickle.load(handle)
            for key in self.def_data_names:
                self.fft_data[key] = pd.read_pickle(str('fft_' + key + '.pickle'))
        for key in self.def_data_names:
            self.fft_data[key] = self.fft_data[key].loc[
                (self.fft_data[key]['freq'] <= self.upper_freq_lim)
            ].reset_index(drop=True)
        self._determine_real_el_freq()
        if fault_display:
            self.fault_freq, self.fault_freq_style = self.fault_freq_detection(fault_display)

    def _extended_park_transformation(self):
        data = self.numpy_data
        park_data = {}
        start = time.time()
        for key in self.def_data_names:
            key_time = time.time() - start
            print(key_time)
            park_data[key] = self._park_transformation(data[key])
        return park_data

    def _park_transformation(self, data):
        print('ok')
        start = time.time()
        i_d = np.sqrt(2/3) * data[0] - np.sqrt(1/6) * data[1] - np.sqrt(1/6) * data[2]
        id_time = time.time() - start
        print(id_time)
        i_q = np.sqrt(1/2) * data[1] - np.sqrt(1/1) * data[2]
        iq_time = time.time() - start
        print(iq_time)
        i_p = np.sqrt(i_d**2 + i_q**2)
        ip_time = time.time() - start
        print(ip_time)
        print('---')
        park_data = {
            'id': i_d,
            'iq': i_q,
            'ip': i_p,
        }
        return park_data

    def run_fft(self):
        """Run fft algorithm on data

        Optionally run when creating FaultAnalyzer.
        Additionally calculating base frequency for each data series.
        Rewriting 'fft_data.pickle'.
        """
        data = self.numpy_park_data
        fft_data = {}
        for key in self.def_data_names:
            tmp = data[key]
            tmp_fft = []
            le = self.data_len[key]
            fft_data[key] = {}
            for idx in range(3):
                fft = sp.fft(tmp[idx]) / le
                fft = abs(fft[range(int(le / 2))])
                tmp_fft.append(fft[:])
                key2 = list(self.park_data[key].keys())
                fft_data[key][key2[idx]] = tmp_fft[idx]
            fft_data[key]['freq'] = np.fft.fftfreq(self.data_len[key], d=self.data_sp[key])[range(int(le / 2))]
        return fft_data

    def _determine_real_el_freq(self):
        """

        Calculate base frequency of data within given frequency range of rated electric frequency
        """
        tmp_freq_range = 5
        tmp_freq = self.freq_el['h']
        for key in self.fft_data:
            tmp_data = self.fft_data[key].loc[
                self.fft_data[key]['freq'].between(tmp_freq - tmp_freq_range, tmp_freq + tmp_freq_range)
            ]
            mx = tmp_data.filter(regex='i')[2:].idxmax()
            self.freq_el[key] = tmp_data['freq'].loc[
                int(np.mean(mx))
            ]

    def fault_freq_detection(self, fault_display):
        faults = {
            'bearing': self._fault_bearing,
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
        """Outer raceway fault freq based on Chapter4/Slide20"""
        tmp = {}
        for key in self.def_data_names:
            i = 1
            val = 0
            tmp[key] = []
            while val < self.upper_freq_lim:
                val = self.freq_el[key] * i
                tmp[key].append(val)
                i += 1
            # tmp[key] = [self.freq_el[key] * k for k in range(1, 8)]
        return tmp

    def _fault_bearing(self):
        """Bearing Fault frequencies based on Chapter4/Slide20."""
        K = 0.4
        N_b = 9
        f_c = K * N_b * self.freq_mech
        tmp = {}
        for key in self.def_data_names:
            f_s = self.freq_el[key]
            i = 1
            val = 0
            tmp[key] = []
            while val < self.upper_freq_lim:
                val = f_s + pow(-1, i) * int((i+1)/2) * f_c
                tmp[key].append(val)
                i += 1
                # tmp[key] = [f_s + pow(-1, x) * int((x+1)/2) * f_c for x in range(1, 16)]
        return tmp

    def plot_all_faults(
        self,
        fault_freq_display=None,
        x_limit=100,
    ):
        """Plotting all available fft data in corresponding graph.

        Parameters
        ----------
        fault_freq_display: bool
            Display fault frequencies from 'self.fault_freq'
        """
        data_plot = {}
        data = self.fft_data
        mx = {}
        mi = {}
        fig = {}
        ax1 = {}
        ax2 = {}
        for key in self.def_data_names:
            fig[key], (ax1[key], ax2[key]) = plt.subplots(2, 1)
            data_plot[key] = data[key].copy()
            mx[key] = data_plot[key].iloc[2:, 1:].max().max()
            mi[key] = data_plot[key].min().min()
            data_plot[key].plot(
                ax=ax1[key],
                x=list(data_plot[key].filter(regex='freq'))[0],
                y=list(data_plot[key].filter(regex='i')),
                title=str(key),
            )
            data_plot[key].plot(
                ax=ax2[key],
                x=list(data_plot[key].filter(regex='freq'))[0],
                y=list(data_plot[key].filter(regex='i')),
                logy=True
            )
            ax2[key].set_xlabel('frequencies')
        for ax in ax1, ax2:
            for key in ax.keys():
                # ax[key].legend('upper right'),
                ax[key].grid('on')
                ax[key].set_xlim([-5, x_limit])
                ax[key].set_ylim([None, 1.2 * mx[key]])
                if fault_freq_display:
                    self._display_fault_freq(ax[key])
        plt.show()

    def plot_raw_data(self, data_names=None):
        data_plot = {}
        data = self.data
        mx = {}
        mi = {}
        fig = {}
        ax1 = {}
        ax2 = {}
        for key in data_names:
            fig[key], (ax1[key], ax2[key]) = plt.subplots(2, 1)
            data_plot[key] = data[key].copy()
            mx[key] = data_plot[key].iloc[2:, 1:].max().max()
            mi[key] = data_plot[key].min().min()
            data_plot[key].plot(
                ax=ax1[key],
                x=list(data_plot[key].filter(regex='t'))[0],
                y=list(data_plot[key].filter(regex='i')),
                title=str(key),
            )
            data_plot[key].plot(
                ax=ax2[key],
                x=list(data_plot[key].filter(regex='t'))[0],
                y=list(data_plot[key].filter(regex='i')),
                logy=True
            )
            ax2[key].set_xlabel('frequencies')
        for ax in ax1, ax2:
            for key in ax.keys():
                # ax[key].legend('upper right'),
                ax[key].grid('on')
                # ax[key].set_xlim([-5, x_limit])
                # ax[key].set_ylim([None, 1.2 * mx[key]])
        plt.show()

    def plot_fault_in_one(
        self,
        item='ia',
        data_names=('h', 'f1', 'f2', 'f3'),
        x_limit=100,
        fault_freq_display=False,
    ):
        """Plotting selected phases of all assigned fft data sets into one graph.

        Parameters
        ----------
        item: str
            Phase (current) to be displayed.
        data_names: list of str
            Data series (faults and/or healthy) to be displayed.
        x_limit: float
            Outer right limit for x axis.
        fault_freq_display: bool
            Display fault frequencies from 'self.fault_freq'
        """
        fig, (ax1, ax2) = plt.subplots(2, 1)
        data_plot = {}
        data = self.fft_data
        mx = 0
        mi = 0
        for key in data_names:
            # data_plot[key] = pd.DataFrame({
            #     str('freq_' + key): data[key]['freq'].loc[(data[key]['freq'] <= 500)],
            #     str(key + '_' + item): data[key][item].loc[(data[key]['freq'] <= 500)],
            # })
            data_plot[key] = pd.DataFrame({
                str('freq_' + key): data[key]['freq'],
                str(key + '_' + item): data[key][item],
            })
            mx_tmp = data[key][item][2:].max()
            mi_tmp = data[key][item].min()
            mx = mx_tmp if mx_tmp > mx else mx
            mi = mi_tmp if mi_tmp < mi else mi
        for key in data_plot:
            data_plot[key].plot(
                ax=ax1,
                x=list(data_plot[key].filter(regex='freq'))[0],
                y=list(data_plot[key].filter(regex=item))[0],
                color=self.def_data_names_col[key],
            )
            data_plot[key].plot(
                ax=ax2,
                x=list(data_plot[key].filter(regex='freq'))[0],
                y=list(data_plot[key].filter(regex=item))[0],
                color=self.def_data_names_col[key],
                logy=True
            )
        for ax in ax1, ax2:
            for key in self.def_data_names:
                self._additional_plot_instructions(
                    ax,
                    {key: self.freq_el[key]},
                    custom_col=self.def_data_names_col[key]
                )
            ax.grid('on')
            ax.set_xlim([-5, x_limit])
            ax.set_ylim([None, 1.2 * mx])
            if fault_freq_display:
                self._display_fault_freq(ax)
            ax.legend(loc='upper right')
        ax2.set_xlabel('frequencies')
        plt.show()
        
    def _display_fault_freq(self, ax):
        for item in self.fault_freq.keys():
            for key, val in self.fault_freq[item].items():
                self._additional_plot_instructions(
                    ax,
                    {str(item + '_' + key): val},
                    label_text=True,
                    custom_col=self.def_data_names_col[key],
                    custom_style=self.fault_freq_style[item],
                )

    def _additional_plot_instructions(self, ax, *data, label_text=False, custom_col=None, custom_style=None):
        """Plot related, manages additional plotting instructions."""
        for key, value in data[0].items():
            self._axvlines(
                value,
                ax=ax,
                color='k' if custom_col is None else custom_col,
                linestyle='solid' if custom_style is None else custom_style,
                lw=.5,
                label=key if label_text else None,
            )

    def _axvlines(self, xs, ax=None, **plot_kwargs):
        """Plot related, creates vertical indicator lines."""
        xs = sp.array((xs,) if sp.isscalar(xs) else xs, copy=False)
        lims = ax.get_ylim()
        x_points = sp.repeat(xs[:, None], repeats=3, axis=1).flatten()
        y_points = sp.repeat(sp.array(lims + (sp.nan,))[None, :],
                             repeats=len(xs), axis=0).flatten()
        plot = ax.plot(x_points, y_points, scaley=False, **plot_kwargs)
        return plot



