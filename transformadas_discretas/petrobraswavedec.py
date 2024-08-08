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
