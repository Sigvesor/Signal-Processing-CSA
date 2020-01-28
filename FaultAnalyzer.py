"""
Created on Mon Jan 22 12:14:10 2020

@author: Ernst
"""

import pandas as pd
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import pickle


class FaultAnalyzer:
    def __init__(self, run_fft=False):
        self.fault1 = pd.read_pickle('fault1.pickle')
        self.fault2 = pd.read_pickle('fault2.pickle')
        self.fault3 = pd.read_pickle('fault3.pickle')
        self.healthy = pd.read_pickle('data_healthy.pickle')
        self.def_data_names = ['h', 'f1', 'f2', 'f3']

        # fft related
        self.fft_data = {}

        self.freq_el = {
            key: 48.8 for key in self.def_data_names
        }
        self.freq_mech = 1420/60
        self.fault_vib = {
            'inner': 5.4152,  # inner ring
            'outer': 3.5848,  # outer ring
            'cage': 0.39828,  # cage train
            'roll': 4.7135,  # rolling element
        }

        h_time = self.healthy['t']
        f1_time = self.fault1['t']
        f2_time = self.fault2['t']
        f3_time = self.fault3['t']
        self.time = {
            'h': h_time,
            'f1': f1_time,
            'f2': f2_time,
            'f3': f3_time,
        }
        h_n = h_time.count()
        f1_n = f1_time.count()
        f2_n = f2_time.count()
        f3_n = f3_time.count()
        self.data_n = {
            'h': h_n,
            'f1': f1_n,
            'f2': f2_n,
            'f3': f3_n,
        }
        h_sp = h_time[1]
        f1_sp = f1_time[1]
        f2_sp = f2_time[1]
        f3_sp = f3_time[1]
        self.data_sp = {
            'h': h_sp,
            'f1': f1_sp,
            'f2': f2_sp,
            'f3': f3_sp,
        }
        # print('         sample time  | sample length \n' +
        #       'healthy |   ' + str(h_sp) + '     |  ' + str(h_n) + '\n' +
        #       'fault1  |   ' + str(f1_sp) + '     |  ' + str(f1_n) + '\n' +
        #       'fault2  |   ' + str(f2_sp) + '     |  ' + str(f2_n) + '\n' +
        #       'fault3  |   ' + str(f3_sp) + '     |  ' + str(f3_n))
        if run_fft:
            self.run_fft()
        else:
            with open('fft_data.pickle', 'rb') as handle:
                self.fft_data, self.freq_el = pickle.load(handle)

    def run_fft(self):
        """Run fft algorithm on data


        """
        data = {
            'h': self.healthy.transpose().iloc[1:].to_numpy(),
            'f1': self.fault1.transpose().iloc[1:].to_numpy(),
            'f2': self.fault2.transpose().iloc[1:].to_numpy(),
            'f3': self.fault3.transpose().iloc[1:].to_numpy(),
        }
        for key in self.def_data_names:
            tmp = data[key]
            tmp_fft = []
            le = self.data_n[key]
            for idx in range(3):
                fft = sp.fft(tmp[idx]) / le
                fft = abs(fft[range(int(le / 2))])
                tmp_fft.append(fft[:])
            self.fft_data[key] = pd.DataFrame({
                'freq': np.fft.fftfreq(self.data_n[key], d=self.data_sp[key])[range(int(le / 2))],
                'ia': tmp_fft[0],
                'ib': tmp_fft[1],
                'ic': tmp_fft[2],
            })
            self.fft_data[key] = self.fft_data[key].loc[
                (self.fft_data[key]['freq'] <= 500)
            ].reset_index(drop=True)
        self._determine_real_el_freq()
        with open('fft_data.pickle', 'wb') as handle:
            pickle.dump((self.fft_data, self.freq_el), handle)

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

    def plot_all_faults(self, logy=True):
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
            for item in data_plot[key].keys():
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
                ax[key].legend('upper right'),
                ax[key].grid('on')
                ax[key].set_xlim([-5, 100])
                ax[key].set_ylim([None, 1.2 * mx[key]])
        plt.show()

    def plot_fault_in_one(self, item='ia', data_names=('h', 'f1', 'f2', 'f3'), x_limit=100, fault_vib_label=False):
        """Plotting selected phases of all assigned data sets into one graph.

        Parameters
        ----------
        item: str
            Phase (current) to be displayed.
        data_names: list of str
            Data series (faults and/or healthy) to be displayed.
        x_limit: float
            Outer right limit for x axis.
        fault_vib_label: bool
            Display faulty vibration legend labels
        """
        fig, (ax1, ax2) = plt.subplots(2, 1)
        data_plot = {}
        data = self.fft_data
        mx = 0
        mi = 0
        for key in data_names:
            data_plot[key] = pd.DataFrame({
                str('freq_' + key): data[key]['freq'].loc[(data[key]['freq'] <= 500)],
                str(key + '_' + item): data[key][item].loc[(data[key]['freq'] <= 500)],
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
            )
            data_plot[key].plot(
                ax=ax2,
                x=list(data_plot[key].filter(regex='freq'))[0],
                y=list(data_plot[key].filter(regex=item))[0],
                logy=True
            )
        for ax in ax1, ax2:
            self._additional_plot_instructions(
                ax,
                {
                    **self.freq_el,
                    # **{key: val * self.freq_mech for key, val in self.fault_vib.items()},
                },
                label_text=fault_vib_label,
            )
            ax.legend(loc='upper right')
            ax.grid('on')
            ax.set_xlim([-5, x_limit])
            ax.set_ylim([None, 1.2 * mx])
        ax2.set_xlabel('frequencies')
        plt.show()

    def _additional_plot_instructions(self, ax, *data, label_text=False):
        """Plot related, manages additional plotting instructions."""
        for key, value in data[0].items():
            self._axvlines(
                value,
                ax=ax,
                # color='r',
                linestyle='--',
                lw=.8,
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



