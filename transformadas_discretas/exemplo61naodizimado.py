import matplotlib.pyplot as plt
import pywt
import numpy as np
import csv

array = [1, -2, 0, 3, 1, 2, -1, 2]
coeffs = pywt.swt(array, 'db1', level=3, trim_approx=True)
fig, axs = plt.subplots(5)
fig.suptitle('Exemplo 6.1 teste')
axs[0].plot(array)
axs[1].plot(coeffs[0])
print(coeffs)
axs[1].set_title('Aproximação')
for i in range (3):
    axs[i+2].set_title('Nível '+ str(i))
    axs[i+2].stem(coeffs[i+1])
plt.tight_layout()
plt.show()
