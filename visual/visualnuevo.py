import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



#### CARGAMOS LOS DATOS DE PM
f = lambda x : (x.replace(",","."))
# PM 10
data = pd.read_csv("datosU.csv", sep = ",", encoding = "ISO=8859-1", decimal = ".")
result = data.values.tolist()
arriba = []
abajo = []

#ESCOGEMOS LAS CELDAS
# celdas del 0 al 7
cell1 = 1 # arriba
cell2 = 3 # abajo

#recorro los distintos pasos de tiempo t de la data
for t in result:
    value = t[cell1]
    value2 = t[cell2]
    arriba.append(value) #este se vera arriba
    abajo.append(value2)

largo = len(arriba)
x = np.arange(0, largo)

largo2 = len(abajo)
x2 = np.arange(0, largo2)

print(arriba)

fig, axs = plt.subplots(2, sharex = True)
fig.suptitle('title')
axs[0].plot(x, arriba)
axs[0].set_title('PM 10')
axs[0].set_ylabel('Concentracion')
axs[1].plot(x2, abajo, 'tab:red')
axs[1].set_title('PM 2.5')
axs[1].set_ylabel('Concentracion')
plt.xlabel('tiempo t')
plt.show()