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
coeffs = pywt.wavedec(retornos_final, 'db1', level=11)
fig, axs = plt.subplots(8)
fig.suptitle('retornos diarios açao petrobras')
axs[0].plot(retornos_final)
axs[1].plot(coeffs[0])
axs[1].set_title('Aproximação')
for i in range (6):
    axs[i+2].set_title('Nível '+ str(i))
    axs[i+2].stem(coeffs[i+1])

plt.show()