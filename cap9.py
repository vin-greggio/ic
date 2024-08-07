
import scipy
import numpy
import math
import pywt
import scipy.signal
from spectrum import data_cosine, Periodogram
import matplotlib.pyplot as plt
from spectrum import periodogram

#def faz_integr_pi(l, jmax)
def faz_H(T, k):
    ans = 0
    for i in range(T):
        ans = ans + numpy.power((1 - numpy.cos(2*numpy.pi*(i/T))), k)
    return ans
    
def closest_value(my_list, target_value):
    return min(my_list, key=lambda x: abs(target_value - x))

def periodizar(psi,x):
    #somatorio de 0 a len(x)
    res = []
    maximo = len(x)
    soma = 0
    dominio = []
    for i in range(maximo):
        soma += 1


data = data_cosine(N=1024, sampling=1024, freq=132.5)
coeffs = pywt.wavedec(data, 'dgau8', mode='periodic', level=3)
plt.plot(coeffs)

wavelet = pywt.Wavelet('sym10')
[phi, psi, x] = wavelet.wavefun(6)
print(len(x), x)
# periodização, dilação e translação.

l = 3
j_max = 5 #the empirical coefficients being zeroed for j>5
a = faz_H(1024, 2)
b = faz_H(1024, 4)
integr_pi = 0

estimador_var = 2*numpy.pi*b*integr_pi/(a*a)
