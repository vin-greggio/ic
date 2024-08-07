import matplotlib.pyplot as plt
import pywt
import numpy as np
import csv

with open('PBR Historical Data (1).csv', encoding = "utf_8") as csvfile:
    reader = csv.reader(csvfile)
    fechamento = []
    abertura = []
    for row in reader:
        fechamento.append(row[1])
        abertura.append(row[2])

retornos = []
for i in range (1,len(abertura)):
    a = float(fechamento[i]) - float(abertura[i])
    retornos.append(a)
retornos_final = retornos[500:]
print(pywt.dwt_max_level(2048,'db8'))
coeffs = pywt.wavedec(retornos_final, 'db4', level=8)
fig, axs = plt.subplots(2,4)
fig.suptitle('retornos diarios açao petrobras')
axs[0,0].plot(retornos_final)
axs[0,1].stem(coeffs[0])
axs[0,1].set_title('Aproximação')
axs[0,2].stem(coeffs[1])
axs[0,3].stem(coeffs[2])
axs[1,0].stem(coeffs[3])
axs[1,1].stem(coeffs[4])
axs[1,2].stem(coeffs[5])
axs[1,3].stem(coeffs[6])


plt.show()
