"""
Created on Mon Jan 22 12:14:10 2020

@author: Ernst
"""

import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt


class FaultAnalyzer:
    def __init__(self):
        self.fault1 = pd.read_pickle('fault1.pickle')
        self.fault2 = pd.read_pickle('fault2.pickle')
        self.fault3 = pd.read_pickle('fault3.pickle')
        self.healthy = pd.read_pickle('data_healthy.pickle')
        self.data_names = ['h', 'f1', 'f2', 'f3']
        self.data_fft = {}

        h_time = self.healthy['t']
        f1_time = self.fault1['t']
        f2_time = self.fault2['t']
        f3_time = self.fault3['t']
        self.time = {
            'h': h_time,
            'f1': f1_time,
            'f2': f1_time,
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
              'fault2   |   ' + str(f2_sp) + '     |  ' + str(f2_n) + '\n' +
              'fault3  |   ' + str(f3_sp) + '     |  ' + str(f3_n))

        self.run_fft()

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
                'freq': sp.arange(int(le / 2)) / le / self.data_sp[item],
                'ia': tmp_fft[0],
                'ib': tmp_fft[1],
                'ic': tmp_fft[2],
            })
            self.data_fft[item] = self.data_fft[item].loc[
                (self.data_fft[item]['freq'] <= 500)
            ].reset_index(drop=True)

    def plot_all_faults(self):
        plt.figure()
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
            'freq': self.data_fft['h']['freq'].loc[(self.data_fft['h']['freq'] <= 500)],
            str('h_' + item): self.data_fft['h'][item].loc[(self.data_fft['h']['freq'] <= 500)],
            str('f1_' + item): self.data_fft['f1'][item].loc[(self.data_fft['f1']['freq'] <= 500)],
            str('f2_' + item): self.data_fft['f2'][item].loc[(self.data_fft['f2']['freq'] <= 500)],
            str('f3_' + item): self.data_fft['f3'][item].loc[(self.data_fft['f3']['freq'] <= 500)],
        })
        mx = data_plot.iloc[2:, 1:4].max().max()
        data_plot.plot(
            ax=ax,
            x='freq',
            y=data_plot.keys()[1:],
        )
        ax.set_xlabel('frequencies')
        ax.legend(loc='upper right')
        ax.grid('on')
        plt.xlim([-5, 100])
        plt.ylim([-.005, 1.2 * mx])
        plt.show()
