"""
Created on Mon Jan 22 13:34:43 2020

@author: Ernst
"""

from FaultAnalyzer import FaultAnalyzer

tmp = FaultAnalyzer(run_fft=False)
# tmp.plot_fault_in_one(
#     item='ic',
#     x_limit=100,
#     fault_vib_label=False,)
    # data_names=['h', 'f1'],
# )
tmp.plot_all_faults()