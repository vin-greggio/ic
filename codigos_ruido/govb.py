#Aqui, foram desenvolvidos métodos de seleção de threshold usando, também, como base a Wavelet Toolbox do MATLAB.
#Neste script, é utilizada a transformada ondaleta discreta e redução de ruído num conjunto de dados que, apesar de possuir
#dados não-equispaçados no tempo, estão numa distribuição aproximadamente uniforme.

import matplotlib.pyplot as plt
import pywt
import numpy as np
import csv
import math

def trivial(var, coeffs, a):
    return abs(a)

def sure_shrink(var, signal, a):
    m = signal.shape[0]
    sorted_signal = np.sort(np.abs(signal))**2
    c = np.linspace(m-1, 0, m)
    s = np.cumsum(sorted_signal) + c * sorted_signal
    risk = (m - (2.0 * np.arange(m)) + s) / m
    ibest = np.argmin(risk)
    thr = np.sqrt(sorted_signal[ibest])
    return thr

def visu_shrink(var, coeffs, a):
    N = len(coeffs)
    thre = math.sqrt(var) * math.sqrt(2 * math.log(N))
    print('vvvv')
    return thre

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


def tsd(th, data, method, mode='soft', wavelets_name='sym8', level=5):
    '''

    :param data: signal
    :param method: {'visushrink', 'sureshrink', 'heursure', 'minmax'}, 'sureshrink' as default
    :param mode: {'soft', 'hard', 'garotte', 'greater', 'less'}, 'soft' as default
    :param wavelets_name: wavelets name in PyWavelets, 'sym8' as default
    :param level: deconstruct level, 5 as default
    :return: processed data
    '''
    methods_dict = {'visushrink': visu_shrink, 'sureshrink': sure_shrink, 'minmax': mini_max, 'trivial': trivial}
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
