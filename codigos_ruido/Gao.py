#Aqui, foram desenvolvidos métodos de seleção de threshold usando, também, como base a Wavelet Toolbox do MATLAB.
#Implementou-se o algoritmo de Gao para estimar o espectro de Fourier de uma função cosseno com ruído.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import random
import statsmodels.api as sm
import csv
import statistics
import pywt
from scipy import signal
from spectrum import data_cosine, pdaniell

def trivial(var, coeffs, a, idx):
    return abs(a)

def gao(var,coeffs,a,idx):
    if idx==1:
        alp = 1.29
    elif idx==2:
        alp = 1.09
    elif idx==3:
        alp = 0.92
    elif idx==4:
        alp = 0.65
    elif idx==5:
        alp = 0.77
    elif idx==6:
        alp = 0.54
    elif idx==7:
        alp = 0.46
    elif idx==8:
        alp = 0.39
    elif idx==9:
        alp = 0.32
    elif idx==10:
        alp = 0.27
    if idx in [1,2,3,4]:
        return math.log(len(data)/2)*alp
    else:
        return math.sqrt(math.log(len(data)))*(math.pi/math.sqrt(3))
    

def sure_shrink(var, signal, a):
    m = signal.shape[0]
    sorted_signal = np.sort(np.abs(signal))**2
    c = np.linspace(m-1, 0, m)
    s = np.cumsum(sorted_signal) + c * sorted_signal
    risk = (m - (2.0 * np.arange(m)) + s) / m
    ibest = np.argmin(risk)
    thr = np.sqrt(sorted_signal[ibest])
    return thr


def visu_shrink(var, coeffs, a, idx):
    N = len(coeffs)
    thre = math.sqrt(var) * math.sqrt(2 * math.log(N))
    return thre

def mini_max(var, coeffs, a, idx):
    N = len(coeffs)
    if N > 32:
        return math.sqrt(var) * (0.3936 + 0.1829 * math.log2(N))
    else:
        return 0
    
def get_var(cD):
    coeffs = cD
    abs_coeffs = []
    for coeff in coeffs:
        abs_coeffs.append(math.fabs(coeff))
    abs_coeffs.sort()
    pos = math.ceil(len(abs_coeffs) / 2)
    var = abs_coeffs[pos] / 0.6545
    return var

def tsd(th, data, method, mode='soft', wavelets_name='sym8', level=5):

    methods_dict = {'visushrink': visu_shrink, 'sureshrink': sure_shrink, 'minmax': mini_max, 'trivial': trivial, 'gao': gao}
    wave = pywt.Wavelet(wavelets_name)

    data_ = data[:]

    (cA, cD) = pywt.dwt(data=data_, wavelet=wave)
    var = get_var(cD)

    coeffs = pywt.wavedec(data=data, wavelet=wavelets_name, level=level)

    for idx, coeff in enumerate(coeffs):
        if idx == 0:
            continue
        thre = methods_dict[method](var, coeff, th, idx)
        coeffs[idx] = pywt.threshold(coeffs[idx], thre, mode=mode)

    thresholded_data = pywt.waverec(coeffs, wavelet=wavelets_name)

    return thresholded_data

def cross_validation(Dados, Wavelet, Modo):
    #Função cedida por Vitor Ribas Perrone
    Even = []
    Odd = []
    for i in range(len(Dados)):
        if i%2 == 0:
            Odd.append(Dados[i])
        else:
            Even.append(Dados[i])
    yOdd = []
    for i in range(int((len(Dados)/2))-1):
        yOdd.append((Odd[i]+Odd[i+1])/2)
    yOdd.append((Odd[i]+Odd[-1])/2)
    yEven = []
    for i in range(int((len(Dados)/2))-1):
        yEven.append((Even[i]+Even[-1])/2)
    yEven.append((Even[i]+Even[-1])/2)
    EvenCoeffs = pywt.wavedec(Even, Wavelet, level = pywt.dwt_max_level(len(Even), Wavelet))
    OddCoeffs = pywt.wavedec(Odd, Wavelet, level = pywt.dwt_max_level(len(Odd), Wavelet))
    Coeffs = []
    for Coeff in EvenCoeffs:
        for coeff in Coeff:
            Coeffs.append(coeff)
    for Coeff in OddCoeffs:
        for coeff in Coeff:
            Coeffs.append(coeff)
    t = 0
    Min = []
    for coeff in Coeffs:
        M = 0
        fEven = tsd(np.abs(coeff), Even, 'trivial', Modo, Wavelet, level = 4)
        fOdd = tsd(np.abs(coeff), Odd, 'trivial', Modo, Wavelet, level = 4)
        for j in range(len(Odd)):
            M = M + (fEven[j]-yOdd[j])** 2 + (fOdd[j]-yEven[j]) **2
        if Min == [] or M < Min[-1]:
            Min.append(M)
            t = np.abs(coeff)
    return np.sqrt(1/(1-np.log(2)/np.log(len(data))))*t


data = data_cosine(N=1024, A=0.05)
f, Pxx_den = signal.periodogram(data, fs=1)
print(f,Pxx_den)



'''
rng = np.random.default_rng()
fs = 10e3
N = 1e5
amp = 2*np.sqrt(2)
freq = 1234.0
noise_power = 0.001 * fs / 2
time = np.arange(N) / fs
data = amp*np.sin(2*np.pi*freq*time)
data += rng.normal(scale=np.sqrt(noise_power), size=time.shape)

f, Pxx_den = signal.periodogram(data, fs)
'''
for i in range(len(Pxx_den)):
    if i == 0:
        Pxx_den[i] = -6
        continue
    Pxx_den[i] = math.log(Pxx_den[i])



denoise1 = tsd(0,Pxx_den, method='gao', mode='soft', wavelets_name='db1', level = pywt.dwt_max_level(len(f),'db1'))

print(pywt.dwt_max_level(len(f),'db1'))


fig, ax = plt.subplots(1,2)
ax[0].plot(Pxx_den)
ax[1].plot(np.exp(denoise1))
plt.show()
