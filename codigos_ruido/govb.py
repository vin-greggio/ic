import matplotlib.pyplot as plt
import pywt
import numpy as np
import csv
import math

def trivial(var, coeffs, a):
    return abs(a)

def ogden(var, coeffs, a):
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
    

def sure_shrink(var, coeffs, a):
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

def visu_shrink(var, coeffs, a):
    N = len(coeffs)
    thre = math.sqrt(var) * math.sqrt(2 * math.log(N))
    print('vvvv')
    return thre

def heur_sure(var, coeffs, a):
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

def mini_max(var, coeffs, a):
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
    methods_dict = {'visushrink': visu_shrink, 'sureshrink': sure_shrink, 'heursure': heur_sure, 'minmax': mini_max, 'ogden': ogden, 'trivial': trivial}
    wave = pywt.Wavelet(wavelets_name)

    data_ = data[:]

    (cA, cD) = pywt.dwt(data=data_, wavelet=wave)
    var = get_var(cD)

    coeffs = pywt.wavedec(data=data, wavelet=wavelets_name, level=level)

    for idx, coeff in enumerate(coeffs):
        if idx == 0:
            continue
        thre = methods_dict[method](var, coeff, th)
        coeffs[idx] = pywt.threshold(coeffs[idx], thre, mode=mode)

    thresholded_data = pywt.waverec(coeffs, wavelet=wavelets_name)

    return thresholded_data

with open('antigos/GOVB.csv', encoding = "utf_8") as csvfile:
    reader = csv.reader(csvfile)
    fechamento = []
    abertura = []
    for row in reader:
        if row[4] == 'null':
            row[4] = 0
        if row[1] == 'null':
            row[1] = 0
        fechamento.append(row[4])
        abertura.append(row[1])

retornos = []
for i in range (2,len(abertura)):
    a = (float(fechamento[i])-float(fechamento[i-1]))/float(fechamento[i])
    retornos.append(math.sqrt(a*a))

final = tsd(0, retornos, 'sureshrink', 'soft', 'db4', level = pywt.dwt_max_level(2048, 'db4'))

plt.plot(final)
plt.show()
