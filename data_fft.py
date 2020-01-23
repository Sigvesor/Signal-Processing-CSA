"""
Created on Mon Jan 22 12:14:10 2020

@author: ernst
"""

data = {
    'h': healthy.transpose().iloc[1:].to_numpy(),
    'f1': fault1.transpose().iloc[1:].to_numpy(),
    'f2': fault2.transpose().iloc[1:].to_numpy(),
    'f3': fault3.transpose().iloc[1:].to_numpy(),
}
data_fft = {}
data_freq = {}
for item in data_names:
    tmp = data[item]
    tmp_fft = []
    le = data_n[item]
    for idx in range(3):
        fft = sp.fft(tmp[idx]) / le
        fft = abs(fft[range(int(le / 2))])
        tmp_fft.append(fft[:])
    data_fft[item] = pd.DataFrame({
        'freq': sp.arange(int(le / 2)) / le / data_sp[item],
        'ia': tmp_fft[0],
        'ib': tmp_fft[1],
        'ic': tmp_fft[2],
    })
    data_fft[item] = data_fft[item].loc[
        (data_fft[item]['freq'] <= 500)
    ].reset_index(drop=True)