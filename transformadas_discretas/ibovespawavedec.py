import matplotlib.pyplot as plt
import pywt
import numpy as np
import csv
from statistics import mean

with open('^BVSP.csv', encoding = "utf_8") as csvfile:
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
for i in range (1,len(abertura)):
    a = float(fechamento[i]) - float(abertura[i])
    retornos.append(a)
retornos_final = retornos[2125:]
coeffs = pywt.wavedec(retornos_final,'db4')
(cA, cD) = pywt.dwt(retornos_final,'db4')
print(pywt.dwt_max_level(2048,'db4'))
fig, axs = plt.subplots(3,3)
fig.suptitle('Decomposição em coeficientes de ondaletas dos retornos diários Ibovespa, ondaleta db4 ')

axs[0,0].plot(retornos_final)
axs[0,1].stem(coeffs[1],markerfmt=" ")
axs[0,1].set_title('Nível 1')
axs[0,2].stem(coeffs[2],markerfmt=" ")
axs[0,2].set_title('Nível 2')
axs[1,0].stem(coeffs[3],markerfmt=" ")
axs[1,0].set_title('Nível 3')
axs[1,1].stem(coeffs[4],markerfmt=" ")
axs[1,1].set_title('Nível 4')
axs[1,2].stem(coeffs[5],markerfmt=" ")
axs[1,2].set_title('Nível 5')
axs[2,0].stem(coeffs[6],markerfmt=" ")
axs[2,0].set_title('Nível 6')
axs[2,1].stem(coeffs[7],markerfmt=" ")
axs[2,1].set_title('Nível 7')
axs[2,2].stem(coeffs[8],markerfmt=" ")
axs[2,2].set_title('Nível 8')


plt.tight_layout()
plt.show()
