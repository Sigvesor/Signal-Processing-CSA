import numpy as np
from numpy import pi
from matplotlib import pyplot as plt
from scipy.fftpack import fft

N = 300000  # Samples
fs = 100000  # Sample frequency (Hz)
Ts = 1.0 / float(fs)  # Period
t = np.linspace(0.0, N * Ts, N)  # Time

fund_freqs = [[49], [50], [51]]  # Fundamental frequencies
freqs = fund_freqs  # Fundamental + Harmonics

for idx, val in enumerate(fund_freqs):
    for i in range(3, 14):
        if i % 2 == 1:
            freqs[idx].append(val[0] * i)


signals = [np.zeros(len(t)), np.zeros(len(t)), np.zeros(len(t))]
magn = np.ones(len(t))
# magn[int(len(t)/3):int(2*len(t)/3)] = 1.5
# magn[:int(len(t)/2)] = 0

for idx, val in enumerate(freqs):
    for i in range(7):
        if i == 0:
            signals[idx] += np.sin(2 * pi * val[i] * t)
        else:
            signals[idx] += 0.15 * np.sin(2 * pi * val[i] * t)
            # signals[idx] += magn * 0.15*np.sin(2*pi*val[i]*t)

if 0:
    plt.figure(figsize=(18, 16))

    plt.subplot(3, 1, 1)
    plt.plot(t, signals[0])

    plt.subplot(3, 1, 2)
    plt.plot(t, signals[1])

    plt.subplot(3, 1, 3)
    plt.plot(t, signals[2])
    plt.show()

if 0:
    plt.figure(figsize=(18, 16))
    a = 0.15*np.sin(2*pi*650*t)
    print(max(a))
    plt.plot(a)
    plt.xlim([0, 10])
    plt.show()

if 1:
    fft_1 = 2 * (abs(fft(signals[0] / len(t))))[:int(len(t) / 2)]
    fft_2 = 2 * (abs(fft(signals[1] / len(t))))[:int(len(t) / 2)]
    fft_3 = 2 * (abs(fft(signals[2] / len(t))))[:int(len(t) / 2)]

    def thd(abs_data):
        sq_sum=0.0
        for r in range(abs(fft_1)):
            sq_sum = sq_sum + (abs(fft_1[r])**2)

        sq_harmonics = sq_sum - (max(abs(fft_1)))**2.0
        thd = 100*sq_harmonics**0.5 / max(abs(fft_1))

        return thd

    print("Total Harmonic Distortion(in percent):")
    print(thd(signals[0][1:int(len(signals[0])/2)]))

    #plt.figure(figsize=(18, 16))
    #plt.plot(fft_1, label="1")
    #plt.plot(fft_2, label="2")
    #plt.plot(fft_3, label="3")
    #plt.xlim([0, 2100])
    #plt.legend()
    #plt.show()



