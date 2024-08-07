import matplotlib.pyplot as plt
import pywt
import numpy as np
import csv

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
print(pywt.dwt_max_level(2048,'db1'))
coeffs = pywt.swt(retornos_final, 'db1', trim_approx=True)


fig, axs = plt.subplots(8)
fig.suptitle('Decomposição em coeficientes de ondaletas dos retornos diários Ibovespa')
axs[0].stem(retornos_final)
axs[1].stem(coeffs[0])
axs[1].set_title('Aproximação')
for i in range (6):
    axs[i+2].set_title('Detalhes - nível '+ str(i))
    axs[i+2].stem(coeffs[i+1])
plt.tight_layout()
plt.show()