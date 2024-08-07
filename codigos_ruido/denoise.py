import matplotlib.pyplot as plt
import pywt
import numpy as np
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

def tsd(th, data, metodo, mode='soft', wavelet='sym8', level=5):
    ondaleta = pywt.Wavelet(wavelet)
    data_ = data[:]

    (cA, cD) = pywt.dwt(data=data_, wavelet=ondaleta)
    var = get_var(cD)
    coeffs = pywt.wavedec(data=data, wavelet=ondaleta, level=level)

    for idx, coeff in enumerate(coeffs):
        if idx == 0:
            continue
        #thre = methods_dict[method](var, coeff, th)
        elif metodo == 'ss':
            thre = sure_shrink(0,coeff,0)
        elif metodo == 'vs':
            thre = visu_shrink(var,coeff,0)
        elif metodo == 'mm':
            thre = mini_max(var,coeff,0)
        elif metodo == 'trivial':
            thre = trivial(var,coeff,th)
        coeffs[idx] = pywt.threshold(coeffs[idx], thre, mode=mode)

    thresholded_data = pywt.waverec(coeffs, wavelet=ondaleta)

    return thresholded_data

def cross_validation(Dados, Wavelet, Modo,n):
    #função desenvolvida em conjunto com Vitor Ribas Perrone
    Even = []
    Odd = []
    for i in range(n):
        if i%2 == 0:
            Odd.append(Dados[i])
        elif 1%2==1:
            Even.append(Dados[i])
    yOdd = []
    yEven = []
    for i in range(int((n/2))-1):
        yOdd.append((Odd[i]+Odd[i+1])/2)
    yOdd.append((Odd[i]+Odd[-1])/2)
    for i in range(int((len(Dados)/2))-1):
        yEven.append((Even[i]+Even[-1])/2)
    yEven.append((Even[i]+Even[-1])/2)
    EvenCoeffs = pywt.wavedec(Even, Wavelet, level = pywt.dwt_max_level(len(Even), Wavelet))
    OddCoeffs = pywt.wavedec(Odd, Wavelet, level = pywt.dwt_max_level(len(Odd), Wavelet))
    Coeffs = []
    for i in range(len(EvenCoeffs)):
        for j in EvenCoeffs[i]:
            Coeffs.append(j)
    for Coeff in OddCoeffs:
        for coeff in Coeff:
            Coeffs.append(j)
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
noise = np.random.normal(0, 0.12, 2048)
for i in range(2048):
    noisy_signal.append(signal[i] + noise[i])

tres = cross_validation(noisy_signal, 'sym6', 'soft', len(noisy_signal))
denoised_signal_cv = tsd(tres, noisy_signal, metodo='trivial', mode='soft', wavelet='sym6', level = 6)
denoised_signal_vs = tsd(4, noisy_signal, metodo='vs', mode='soft', wavelet='sym6', level = 6)
denoised_signal_ss = tsd(4, noisy_signal, metodo='ss', mode='soft', wavelet='sym6', level = 6)
fig, ax = plt.subplots(3,2)
denoised_signal_mm = tsd(4, noisy_signal, metodo='mm', mode='soft', wavelet='sym6', level = 6)
#ax[1,1].plot(denoised_signal_ogd)
ax[0,1].plot(denoised_signal_ss)
ax[0,0].plot(denoised_signal_vs)
ax[1,0].plot(denoised_signal_cv)
ax[2,0].plot(denoised_signal_mm)
ax[1,0].title.set_text('cross validation')
ax[0,0].title.set_text('visushrink')
ax[0,1].title.set_text('sureshrink')
ax[2,0].title.set_text('minimax')
ax[2,1].plot(noisy_signal)
ax[2,1].title.set_text('Sinal com ruídos')
plt.tight_layout()
plt.show()
