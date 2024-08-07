import matplotlib.pyplot as plt
import pywt
import numpy as np
import csv
from statistics import mean

values = []
with open('andromeda.csv', encoding = "utf_8") as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        if row[1] != 'Magnitude':
            values.append(float(row[1]))

print(values)
print(len(values))
values = values[32:]
coeffs = pywt.swt(values, 'db1', level=8, trim_approx=True)
fig, axs = plt.subplots(8)
fig.suptitle('ru andromeda')
axs[0].plot(values)
axs[1].plot(coeffs[0])
axs[1].set_title('Aproximação')
for i in range (6):
    axs[i+2].set_title('Nível '+ str(i))
    axs[i+2].stem(coeffs[i+1])
plt.tight_layout()
plt.show()
