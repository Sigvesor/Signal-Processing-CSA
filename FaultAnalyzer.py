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
        self.data_names = ['h', 'f1', 'f2', 'f3']
        self.data_fft = {}
        self.vib = 48.8

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
        print('         sample time  | sample length \n' +
              'healthy |   ' + str(h_sp) + '     |  ' + str(h_n) + '\n' +
              'fault1  |   ' + str(f1_sp) + '     |  ' + str(f1_n) + '\n' +
              'fault2  |   ' + str(f2_sp) + '     |  ' + str(f2_n) + '\n' +
              'fault3  |   ' + str(f3_sp) + '     |  ' + str(f3_n))
        if run_fft:
            self.run_fft()
        else:
            with open('data_fft.pickle', 'rb') as handle:
                self.data_fft = pickle.load(handle)

    def run_fft(self):
        data = {
            'h': self.healthy.transpose().iloc[1:].to_numpy(),
            'f1': self.fault1.transpose().iloc[1:].to_numpy(),
            'f2': self.fault2.transpose().iloc[1:].to_numpy(),
            'f3': self.fault3.transpose().iloc[1:].to_numpy(),
        }
        for item in self.data_names:
            tmp = data[item]
            tmp_fft = []
            le = self.data_n[item]
            for idx in range(3):
                fft = sp.fft(tmp[idx]) / le
                fft = abs(fft[range(int(le / 2))])
                tmp_fft.append(fft[:])
            self.data_fft[item] = pd.DataFrame({
                # 'freq': sp.arange(int(le / 2)) / le / self.data_sp[item],
                'freq': np.fft.fftfreq(self.data_n[item], d=self.data_sp[item])[range(int(le / 2))],
                'ia': tmp_fft[0],
                'ib': tmp_fft[1],
                'ic': tmp_fft[2],
            })
            self.data_fft[item] = self.data_fft[item].loc[
                (self.data_fft[item]['freq'] <= 500)
            ].reset_index(drop=True)
        with open('data_fft.pickle', 'wb') as handle:
            pickle.dump(self.data_fft, handle)

    def plot_all_faults(self):
        # plt.figure()
        ax = {}
        data_plot = {}
        mx = {}
        for item in self.data_names:
            plt.figure()
            ax[item] = plt.axes()
            data_plot[item] = self.data_fft[item]
            mx[item] = self.data_fft[item].iloc[2:, 1:4].max().max()
            data_plot[item].plot(
                ax=ax[item],
                x='freq',
                y=['ia', 'ib', 'ic'],
                title=str(item)
            )
            ax[item].set_xlabel('frequencies')
            ax[item].legend()
            ax[item].grid('on')
            plt.xlim([-5, 100])
            plt.ylim([-.005, 1.1 * mx[item]])
        plt.show()

    def plot_fault_in_one(self, item='ia'):
        plt.figure()
        ax = plt.axes()
        data_plot = pd.DataFrame({
            # 'freq': self.data_fft['h']['freq'].loc[(self.data_fft['h']['freq'] <= 500)],
            'freq_h': self.data_fft['h']['freq'].loc[(self.data_fft['h']['freq'] <= 500)],
            'freq_f1': self.data_fft['f1']['freq'].loc[(self.data_fft['f1']['freq'] <= 500)],
            'freq_f2': self.data_fft['f2']['freq'].loc[(self.data_fft['f2']['freq'] <= 500)],
            'freq_f3': self.data_fft['f3']['freq'].loc[(self.data_fft['f3']['freq'] <= 500)],
            str('h_' + item): self.data_fft['h'][item].loc[(self.data_fft['h']['freq'] <= 500)],
            str('f1_' + item): self.data_fft['f1'][item].loc[(self.data_fft['f1']['freq'] <= 500)],
            str('f2_' + item): self.data_fft['f2'][item].loc[(self.data_fft['f2']['freq'] <= 500)],
            str('f3_' + item): self.data_fft['f3'][item].loc[(self.data_fft['f3']['freq'] <= 500)],
        })
        mx = data_plot.iloc[2:, 4:].max().max()
        for name in self.data_names:
            data_plot.plot(
                ax=ax,
                # x=['freq_h', 'freq_f1', 'freq_f2', 'freq_f3'],
                x=str('freq_' + name),
                y=str(name + '_' + item),
            )
            # plt.hold(True)
        ax.set_xlabel('frequencies')
        ax.legend(loc='upper right')
        ax.grid('on')
        self._additional_plot_instructions(ax, self.vib)
        plt.xlim([-5, 100])
        plt.ylim([-.005, 1.2 * mx])
        plt.show()

    def _additional_plot_instructions(self, ax, *data):
        """Plot related, manages additional plotting instructions."""
        self._axvlines(
            data[0],
            ax=ax,
            color='r',
            linestyle='--',
            lw=.5,
            # label=label_text,
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
