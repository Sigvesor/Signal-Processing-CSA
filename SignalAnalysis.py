import numpy as np
from numpy import pi

import matplotlib as mp
import matplotlib.pyplot as plt
"""Only needed for latex useable plots"""
# size = 30
# pgf_with_latex = {                      # setup matplotlib to use latex for output
#     "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
#     "text.usetex": True,                # use LaTeX to write all text
#     "font.family": 'sans',
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

from scipy.fftpack import fft, fftfreq
import scipy.signal as sg
from scipy.stats import hmean

class SignalAnalysis():
    def __init__(self, fidx=1):
        self.harmonics = np.arange(1, 14, 2) # harmonic range 1,3,5,7,9,11,13
        self.fund_freqs = [49, 50, 51] # Fundamental frequencie
        self.freqs = []
        self.fidx = fidx
        self.fs = None
        self.ts = None
        self.t_end = 1
        self.exp = ['exp1', 'exp2', 'exp3']
        self.t, self.signals = self._create_data()
        self.output = {}
        self.input = {}
        self.system = self._run_simulation()
        
    def _create_data(self):
        """Create input data"""
        n_harmonics = len(self.harmonics)
        for idx, freq in enumerate(self.fund_freqs):
            self.freqs.append([val * freq for val in self.harmonics])
        self.fs = 10*max(max(self.freqs)) # sample frequency (Hz)
        self.ts = float(self.fs)**-1 # sample time
        T = [x ** -1 for x in self.fund_freqs]        
        t = {}
        signals = {}
        for idx1, val1 in enumerate(self.exp):
            if idx1 == 1:
                t_end = self.t_end * 3
            elif idx1 == 2: 
                t_end = self.t_end * 2
            else:
                t_end = self.t_end
            t[val1] = np.arange(0.0, t_end, self.ts)
            signals[val1] = [np.zeros(len(t[val1])) for i in T]  
        
        for idx1, val1 in enumerate(self.exp):
            l = len(t[val1])
            magn = np.ones(l)
            # signals[val1] = []
            if idx1 == 1:
                magn[int(l/3):int(2*l/3)] = 1.5
            elif idx1 == 2:
                magn[int(l/2):] = 0
            for idx2, val2 in enumerate(self.freqs):
                for i in range(n_harmonics):
                    if i == 0:
                        signals[val1][idx2] += magn * np.sin(2*pi*val2[i]*t[val1])
                    else:
                        signals[val1][idx2] += magn * 0.15*np.sin(2*pi*val2[i]*t[val1])
        return t, signals


    def _run_simulation(self, filter='bandpass', fidx=-1, R=1, L=.1):
        """Run simulation according to input"""
        fidx = fidx if fidx+1 else self.fidx
        system = getattr(self, '_' + filter)(R=R, L=L)
        for idx, val in enumerate(self.exp):
            sin = self.signals[val][fidx]
            tin = self.t[val]
            tout, sout, x = sg.lsim(system, sin, tin)
            sin_freq, sin_fft = self._run_fft(tin, sin)
            sout_freq, sout_fft = self._run_fft(tout, sout)
            sin_thdf = self._calc_thdf(sin_freq, sin_fft, self.fund_freqs[fidx])
            sout_thdf = self._calc_thdf(sout_freq, sout_fft, self.fund_freqs[fidx])
            sin_thdi = self._calc_thdi(sin_freq, sin_fft)
            sout_thdi = self._calc_thdi(sout_freq, sout_fft)
            sin_thdr = self._calc_thdr(sin_thdf)
            sout_thdr = self._calc_thdr(sout_thdf)
            self.input[val] = {
                't': tin,
                's': sin,
                'freq': sin_freq,
                'fft': sin_fft,
                'thdf': sin_thdf,
                'thdr': sin_thdr,
                'thdi': sin_thdi,
            }
            self.output[val] = {
                't': tout,
                's': sout,
                'freq': sout_freq,
                'fft': sout_fft,
                'thdf': sout_thdf,
                'thdr': sout_thdr,
                'thdi': sout_thdi,
                }
        return system
        
    def _bandpass(self, L=.1, R=1, fidx=-1):
        """Create bandpass filter system"""
        fidx = fidx if fidx+1 else self.fidx
        w_0 = 2*pi*self.fund_freqs[fidx]
        # space state representation
        A = [[0, 1], [-w_0**2, -R/L]]
        B = [[0], [R/L]]
        C = [0, 1]
        D = 0;
        system = sg.lti(A, B, C, D)
        return system
    
    def _run_fft(self, t, signal):
        sig = 2*(abs(fft(signal/len(t))))[:int(len(t)/2)]
        freq = fftfreq(len(signal), d=self.ts)[range(int(len(t)/2))]
        return freq, sig
    
    def _calc_thdf(self, fft_freq, fft_sig, f=None, fidx=-1):
        fidx = fidx if fidx+1 else self.fidx
        f = f if f else self.fund_freqs[fidx]
        f_h = self.fund_freqs[fidx]
        tmp_s = fft_sig[fft_freq % f_h == 0]
        tmp_f = fft_freq[fft_freq % f_h == 0]
        hpeaks = tmp_s[(tmp_f != 0) & (tmp_f != f)]
        thdf = np.sqrt(sum(hpeaks**2))/fft_sig[fft_freq == f] * 100
        return thdf[0] 
    
    def _calc_thdi(self, fft_freq, fft_sig, fidx=-1):
        fidx = fidx if fidx+1 else self.fidx
        f_1 = self.fund_freqs[fidx]
        fpeak = fft_sig[fft_freq == f_1]
        h = [fft_freq[fft_freq == f_1 * hi] for hi in self.harmonics[1:]]
        hpeaks = [fft_sig[fft_freq == f_1 * hi] for hi in self.harmonics[1:]]
        thdi = [val/fpeak * 100 for val in hpeaks]
        thdi = [i[0] for i in thdi]
        return thdi
    
    def _calc_thdr(self, thdf):
        thdr = thdf/np.sqrt(1 + thdf**2)
        return thdr
    
    def plot_signals(self, exp='exp1'):
        plt.figure(figsize=(18, 16))
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(18,16))
        for idx, ax in enumerate([ax1, ax2, ax3]):
            ax.plot(self.t[exp], self.signals[exp][idx])
            ax.set_xlabel('time in s')
            ax.grid('minor')
        plt.show()

    def plot_bode(self, filter='bandpass', fidx=-1):
        fig= plt.figure(figsize=(24, 12))
        gs = mp.gridspec.GridSpec(2, 2, figure=fig)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, 0])
        ax4 = fig.add_subplot(gs[1, 1])
        fidx = fidx if fidx+1 else self.fidx
        for i in range(len(self.fund_freqs)):
            system = getattr(self, '_' + filter)(fidx=i)
            w, mag, phase = sg.bode(system, n=1000)
            for ax in ax1, ax3:
                ax.semilogx(w, mag)
                ax.set_ylabel('magnitude')
                ax.set_xlabel('frequency (rad)') 
                ax.grid(True)   
            for ax in ax2, ax4:
                ax.set_ylabel('phase ($^\circ$)')        
                ax.set_xlabel('frequency (rad)')   
                ax.semilogx(w, phase)
                ax.set_yticks([-90, -45, 0, 45, 90])
                ax.grid(True)   
            for ax in ax3, ax4:
                ax.grid(True, which='both')
                ax.set_xlim([150, 600])
        ax2.legend(['$\omega_0 = 2\pi' + str(f) + '\mathrm{Hz}$' for f in self.fund_freqs], loc='lower left')
        ax3.set_ylim([-50, 5])
        fig.tight_layout()
        fig.savefig('bode.pgf')
        # fig.savefig('../../sigproc_report/figures/bode/bode.pgf')

            
    def plot_simulation(self, fidx=-1, exp='exp1'):
        fidx = fidx if fidx+1 else self.fidx
        
        tin = self.input[exp]['t']
        tout = self.output[exp]['t']
        sin = self.input[exp]['s']
        sout = self.output[exp]['s']
        sin_freq = self.input[exp]['freq']
        sout_freq = self.output[exp]['freq']
        sin_fft = self.input[exp]['fft']
        sout_fft = self.output[exp]['fft']
        sin_thdf = self.input[exp]['thdf']
        sout_thdf = self.output[exp]['thdf']
        
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(18,16))
        ax1.plot(tin, sin)
        ax1.plot(tout, sout)
        ax1.legend(['in', 'out'])
        ax1.title.set_text('h_1 = ' + str(self.fund_freqs[fidx]))
        
        w, mag, phase = sg.bode(self.system)
        ax2.semilogx(w/(2*pi), phase)
        ax2.title.set_text('bode')
                
        ax3.plot(sin_freq, sin_fft)
        ax3.plot(sout_freq, sout_fft)
        ax3.set_xlim([0, 1000])
        ax3.legend([
            'in_' + str(sin_thdf) + '%',
            'out_' + str(sout_thdf) + '%'
            ])
        
    def distorsions(self, exp='exp1'):
        """Print individual distorsions to console"""
        print('input thdi:\n' + str(self.input[exp]['thdi']))
        print('output thdi:\n' + str(self.output[exp]['thdi']))
        
    def print_plot_for_report(self, fidx=-1, exp='exp1'):
        """Print plots for report"""
        fidx = fidx if fidx+1 else self.fidx
        loc = '../../sigproc_report/'
        f1 = self.fund_freqs[fidx]
        tin = self.input[exp]['t']
        tout = self.output[exp]['t']
        sin = self.input[exp]['s']
        sout = self.output[exp]['s']
        sin_freq = self.input[exp]['freq']
        sout_freq = self.output[exp]['freq']
        sin_fft = self.input[exp]['fft']
        sout_fft = self.output[exp]['fft']
        sin_thdf = round(self.input[exp]['thdf'], 3)
        sout_thdf = round(self.output[exp]['thdf'], 3)
        sin_thdi = self.input[exp]['thdi']
        sout_thdi = self.output[exp]['thdi']
        
        fig = plt.figure(figsize=(24,12))
        gs = mp.gridspec.GridSpec(2, 2, figure=fig)
        ax1 = fig.add_subplot(gs[0, :])
        ax2 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[1, 1])
        
        ax1.plot(tin, sin)
        ax1.plot(tout, sout)
        ax1.text(
            .5, 
            0, 
            '$h_1$ = ' + str(self.fund_freqs[fidx]) + 'Hz',
            bbox={
                'facecolor': 'white',
                'boxstyle': 'round',
                'edgecolor': 'white',
                'alpha': 1},
        )
        ax1.text(-max(tin)/10, max(sin), 'a)', fontweight='bold')
        ax1.set_xlabel('time (s)')
        
        l = int(1/ (self.fund_freqs[fidx] * self.ts) * 1.1)
        if exp == 'exp1':
            tin = tin[-l:-1]
            sin = sin[-l:-1]
            tout = tout[-l:-1]
            sout = sout[-l:-1]
        elif exp == 'exp2':
            i = int(len(tin)/3)*2
            tin = tin[i-l:i-1]
            sin = sin[i-l:i-1]
            tout = tout[i-l:i-1]
            sout = sout[i-l:i-1]
            ax3.set_yticks([-1.5, -1, -.5, 0 , .5, 1, 1.5])
        elif exp == 'exp3':
            i = int(len(tin)/2)
            tin = tin[i-l:i-1]
            sin = sin[i-l:i-1]
            tout = tout[i-l:i-1]
            sout = sout[i-l:i-1]
            
        ax3.plot(tin, sin)
        ax3.plot(tout, sout)
        ax3.text(tin[0]-.0035, max(sin), 'c)', fontweight='bold')
        ax3.set_xlabel('time (s)')

        ax2.plot(sin_freq, sin_fft)
        ax2.plot(sout_freq, sout_fft)
        ax2.set_xlim([0, 750])
        ax2.legend([
            'in (THD=' + str(sin_thdf) + '%)',
            'out (THD=' + str(sout_thdf) + '%)'
            ])
        ax2.set_xlabel('frequency (Hz)')
        ax2.text(-85, max(sin_fft), 'b)', fontweight='bold')
        
        for ax in [ax1, ax2, ax3]:
            ax.grid(True)
            ax.set_ylabel('magnitude')
            
        fig.tight_layout()
        fig.savefig(loc + 'figures/' + exp + '/all' + str(f1) + '.pgf')
        print(
            '\SI{' + str(sin_thdf) + '}{\percent}', 
            file=open(loc + 'values/' + exp + '_thd_in' + str(f1) + '.tex', 'w'),
        )
        print(
            '\SI{' + str(sout_thdf) + '}{\percent}', 
            file=open(loc + 'values/' + exp + '_thd_out' + str(f1) + '.tex', 'w'),
        )
        print(
            ' & '.join(['\SI{' + str(round(val, 3)) + '}{\percent}' for val in sin_thdi]), 
            file=open(loc + 'values/' + exp + '_thdi_in' + str(f1) + '.tex', 'w'),
        )
        print(
            ' & '.join(['\SI{' + str(round(val, 3)) + '}{\percent}' for val in sout_thdi]), 
            file=open(loc + 'values/' + exp + '_thdi_out' + str(f1) + '.tex', 'w'),
        )
        
        
    def plot_comparison_for_report(self, fidx=-1, exp='exp1'):
        """Plotting comparison for mismatching w0 and f1"""
        fidx = fidx if fidx+1 else self.fidx
        loc = '../../sigproc_report/'
        fig = plt.figure(figsize=(24,12))
        gs = mp.gridspec.GridSpec(2, 2, figure=fig)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, 0])
        ax4 = fig.add_subplot(gs[1, 1])
        l1 = int(1/ (self.fund_freqs[fidx] * self.ts) * 30)
        l2 = int(1/ (self.fund_freqs[fidx] * self.ts) * 1.1)
        
        ax = [[ax1, ax3], [ax2, ax4]]
        for i in range(2):
            if i == 1:
                self._run_simulation(R=10)
            else:
                self._run_simulation()
            tin = self.input[exp]['t']
            tout = self.output[exp]['t']
            sin = self.input[exp]['s']
            sout = self.output[exp]['s']
            sin_thdf = round(self.input[exp]['thdf'], 3)
            sout_thdf = round(self.output[exp]['thdf'], 3)
            
            ax[i][0].plot(tin[0:l1], sin[0:l1])
            ax[i][0].plot(tout[0:l1], sout[0:l1])
            ax[i][1].plot(tin[-l2:-1], sin[-l2:-1])
            ax[i][1].plot(tout[-l2:-1], sout[-l2:-1])
            ax[i][1].legend([
            'in (THD=' + str(sin_thdf) + '\%)',
            'out (THD=' + str(sout_thdf) + '\%)'],
            loc='upper right',
            fontsize=20
            )
        for ax in ax1, ax2, ax3, ax4:
            ax.grid(True)
        ax3.set_xlabel('time (s)')
        ax1.text(.3, 1.3, 'a)', fontweight='bold')
        ax2.text(.3, 1.3, 'b)', fontweight='bold')
        ax4.set_xlabel('time (s)')
        ax1.set_ylabel('magnitude')
        ax3.set_ylabel('magnitude')          
        fig.tight_layout()
        fig.savefig(loc + 'figures/' + exp + '/comparison.pgf')

if __name__ == "__main__":
    """Show all three 'exp' for all three 'fund freqs' for report"""
    for i in [0, 1, 2]:
        tmp = SignalAnalysis(fidx=i)
        for idx in tmp.exp:
            tmp.print_plot_for_report(exp=idx)
    
    
    # tmp = SignalAnalysis(fidx=2)
    # for idx in tmp.exp:
    """"Plot Simulation results"""
    #     tmp.plot_simulation(exp=idx)
    """plot specific exp for report"""
    # tmp.print_plot_for_report(exp='exp2')
    """plot design comparison"""
    # tmp.plot_comparison_for_report(fidx=2)
    """plot bode plot"""
    tmp.plot_bode()
    """print distorsions"""
    tmp.distorsions()