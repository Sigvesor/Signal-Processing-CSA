"""
Created on Mon Jan 22 13:34:43 2020

@author: Ernst
"""

from FaultAnalyzer import FaultAnalyzer

fault_display = [
    # 'bearing',
    'stf',
]
tmp = FaultAnalyzer(run_park_tr=True, run_fft=True, fault_display=fault_display, upper_freq_lim=1000,)
tmp.plot_fault_in_one(
    item='ip',
    x_limit=300,
    data_names=[
        # 'h',
        'f1',
        'f2',
        'f3',
    ],
    fault_freq_display=True,
)
# tmp.plot_all_faults(
#     fault_freq_display=True,
#     x_limit=500,
#     data_names=[
#         'h',
#         'f1',
#         'f2',
#         'f3',
#     ],
# )
# tmp.plot_raw_data(
#     data_names=['f1']
# )
tmp.plot_faults_comparison()
