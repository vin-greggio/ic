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
        return math.log(len(data))*alp
    else:
        return math.sqrt(math.log(len(data)))*(math.pi/math.sqrt(3))

def ogden(var, coeffs, a, idx):
    coeffs_dumm = []
    for i in range(len(coeffs)):
        coeffs_dumm.append(abs(coeffs[i]))
    #alpha-quantile
    maxim = math.pow(max(coeffs_dumm), 2)
    for i in range(len(coeffs)):
        squared = math.pow(coeffs_dumm[i], 2)
        if (squared>1.6449):
            coeffs_dumm[i] = 0
    thre = max(coeffs_dumm)
    print(thre)
    return thre
    

def sure_shrink(var, coeffs, a, idx):
    N = len(coeffs)
    sqr_coeffs = []
    for coeff in coeffs:
        sqr_coeffs.append(math.pow(coeff, 2))
    sqr_coeffs.sort()
    pos = 0
    r = 0
    for idx, sqr_coeff in enumerate(sqr_coeffs):
        new_r = (N - 2 * (idx + 1) + (N - (idx + 1))*sqr_coeff + sum(sqr_coeffs[0:idx+1]))
        if r == 0 or r > new_r:
            r = new_r
            pos = idx
    thre = math.sqrt(sqr_coeffs[pos])
    return thre

def visu_shrink(var, coeffs, a, idx):
    N = len(coeffs)
    thre = math.sqrt(var) * math.sqrt(2 * math.log(N))
    return thre

def heur_sure(var, coeffs, a, idx):
    N = len(coeffs)
    s = 0
    for coeff in coeffs:
        s += math.pow(coeff, 2)
    theta = (s - N) / N
    miu = math.pow(math.log2(N), 3/2) / math.pow(N, 1/2)
    if theta < miu:
        return visu_shrink(var, coeffs)
    else:
        return min(visu_shrink(var, coeffs), sure_shrink(var, coeffs))

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

def get_baseline(data, wavelets_name='sym8', level=5):
    '''

    :param data: signal
    :param wavelets_name: wavelets name in PyWavelets, 'sym8' as default
    :param level: deconstruct level, 5 as default
    :return: baseline signal
    '''
    wave = pywt.Wavelet(wavelets_name)
    coeffs = pywt.wavedec(data, wave, level=level)
    for i in range(1, len(coeffs)):
        coeffs[i] *= 0
    baseline = pywt.waverec(coeffs, wave)
    return baseline

def tsd(th, data, method, mode='soft', wavelets_name='sym8', level=5):
    '''

    :param data: signal
    :param method: {'visushrink', 'sureshrink', 'heursure', 'minmax'}, 'sureshrink' as default
    :param mode: {'soft', 'hard', 'garotte', 'greater', 'less'}, 'soft' as default
    :param wavelets_name: wavelets name in PyWavelets, 'sym8' as default
    :param level: deconstruct level, 5 as default
    :return: processed data
    '''
    methods_dict = {'visushrink': visu_shrink, 'sureshrink': sure_shrink, 'heursure': heur_sure, 'minmax': mini_max, 'ogden': ogden, 'trivial': trivial, 'gao': gao}
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