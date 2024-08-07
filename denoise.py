import matplotlib.pyplot as plt
import pywt
import numpy as np
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
    print('ssss')
    return thre

def visu_shrink(var, coeffs, a):
    N = len(coeffs)
    thre = math.sqrt(var) * math.sqrt(2 * math.log(N))
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
        fEven = tsd(np.abs(coeff), Even, 'trivial', Modo, Wavelet, level = 6)
        fOdd = tsd(np.abs(coeff), Odd, 'trivial', Modo, Wavelet, level = 6)
        for j in range(len(Odd)):
            M = M + (fEven[j]-yOdd[j])** 2 + (fOdd[j]-yEven[j]) **2
        if Min == [] or M < Min[-1]:
            Min.append(M)
            t = np.abs(coeff)
    return np.sqrt(1/(1-np.log(2)/np.log(len(noisy_signal))))*t
        

base = []
signal = []
contador = 0
for i in range(2048):
    base.append(contador + 1/2048)
    contador += 1/2048
for i in base:
    signal.append(np.sin((2.1*math.pi)/(i+0.05)))

noisy_signal = []
noise = np.random.normal(0, 0.25, 2048)
for i in range(2048):
    noisy_signal.append(signal[i] + noise[i])

tres = cross_validation(noisy_signal, Wavelet='sym6', Modo='soft')
denoised_signal_cv = tsd(tres, noisy_signal, method='trivial', mode='soft', wavelets_name='sym6', level = 6)
denoised_signal_vs = tsd(4, noisy_signal, method='visushrink', mode='soft', wavelets_name='sym6', level = 6)
denoised_signal_ss = tsd(4, noisy_signal, method='minmax', mode='soft', wavelets_name='sym6', level = 6)
fig, ax = plt.subplots(3,2)
denoised_signal_ogd = tsd(4, noisy_signal, method='ogden', mode='soft', wavelets_name='sym6', level = 6)
denoised_signal_mm = tsd(4, noisy_signal, method='minmax', mode='soft', wavelets_name='sym6', level = 6)
ax[1,1].plot(denoised_signal_ogd)
ax[0,1].plot(denoised_signal_ss)
ax[0,0].plot(denoised_signal_vs)
ax[1,0].plot(denoised_signal_cv)
ax[2,0].plot(denoised_signal_mm)
ax[1,0].title.set_text('cross validation')
ax[0,0].title.set_text('visushrink')
ax[1,1].title.set_text('ogden&parzen')
ax[0,1].title.set_text('sureshrink')
ax[2,0].title.set_text('minimax')
ax[2,1].plot(noisy_signal)
ax[2,1].title.set_text('Sinal com ruídos')
plt.tight_layout()
plt.show()