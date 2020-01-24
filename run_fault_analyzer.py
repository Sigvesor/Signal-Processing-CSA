"""
Created on Mon Jan 22 13:34:43 2020

@author: Ernst
"""

from FaultAnalyzer import FaultAnalyzer

tmp = FaultAnalyzer(run_fft=False)
tmp.plot_fault_in_one('ic')
# tmp.plot_all_faults()