"""
Created on Mon Jan 22 13:34:43 2020

@author: Ernst
"""

from FaultAnalyzer import FaultAnalyzer

fault_display = [
    # 'bearing',
    'stf',
]
tmp = FaultAnalyzer(run_fft=False, fault_display=fault_display, upper_freq_lim=5000,)
tmp.plot_fault_in_one(
    item='ic',
    x_limit=1000,
    data_names=[
        'h',
        'f1',
        'f2',
        'f3',
    ],
    fault_freq_display=True,
)
# tmp.plot_all_faults(
#     fault_freq_display=True,
# )
