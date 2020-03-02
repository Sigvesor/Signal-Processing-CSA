# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 21:25:47 2020

@author: Sigve Sørensen
"""

import pandas as pd
from scipy.io import loadmat

def to_pickle():
# Load .mat files
    fault1_mat = loadmat("./raw_data/fault1.mat")
    fault2_mat = loadmat("./raw_data/fault2.mat")
    fault3_mat = loadmat("./raw_data/fault3.mat")
    healthy_mat = loadmat("./raw_data/data_healthy.mat")

    # Convert to pd DataFrames
    # Correct labels
    data1 = fault1_mat['data']
    fault1 = pd.DataFrame(data1)
    fault1[3] = fault1_mat['time']
    fault1 = fault1[[3, 0, 1, 2]]
    fault1.rename(columns={3: 't',
                           0: 'ia',
                           1: 'ib',
                           2: 'ic'}, inplace=True)

    data2 = fault2_mat['time']
    fault2 = pd.DataFrame(data2)
    fault2['ia'] = fault2_mat['ia']
    fault2['ib'] = fault2_mat['ib']
    fault2['ic'] = fault2_mat['ic']
    fault2.rename(columns={0: 't'}, inplace=True)

    data3 = fault3_mat['time']
    fault3 = pd.DataFrame(data3)
    fault3['ia'] = fault3_mat['ia']
    fault3['ib'] = fault3_mat['ib']
    fault3['ic'] = fault3_mat['ic']
    fault3.rename(columns={0: 't'}, inplace=True)

    data4 = healthy_mat['time']
    healthy = pd.DataFrame(data4)
    healthy['ia'] = healthy_mat['ia']
    healthy['ib'] = healthy_mat['ib']
    healthy['ic'] = healthy_mat['ic']
    healthy.rename(columns={0: 't'}, inplace=True)

    # Read to pickle
    fault1.to_pickle('fault1.pickle')
    fault2.to_pickle('fault2.pickle')
    fault3.to_pickle('fault3.pickle')
    healthy.to_pickle('data_healthy.pickle')
