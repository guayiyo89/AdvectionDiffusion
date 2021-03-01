import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



#### CARGAMOS LOS DATOS DE PM
f = lambda x : (x.replace(",","."))
# PM 10
data = pd.read_csv("visual/PM_estimado.csv", sep = ",", encoding = "ISO=8859-1", decimal = ".")
result = data.values.tolist()
pm10 = []

def title(celda):
    titulo = 'Concentraciones'
    if celda == 13:
        titulo = 'Nielol'
    if celda == 7:
        titulo = 'Las Encinas'
    if celda == 24:
        titulo = 'Padre las Casas'
    return titulo

cell = 48
title = title(cell)

for a in result:
    value = a[cell]
    pm10.append(value)


# PM 25
data = pd.read_csv("visual/PM_exacto.csv", sep = ",", encoding = "ISO=8859-1", decimal = ".")
result = data.values.tolist()
pm25 = []

for a in result:
    value = a[cell]
    pm25.append(value)

largo = len(pm10)
x = np.arange(0, largo)

largo2 = len(pm25)
x2 = np.arange(0, largo2)


fig, axs = plt.subplots(2, sharex = True)
fig.suptitle('Evaluacion de la Solucion Manufacturada (Celda 49)')
axs[0].plot(x, pm10)
axs[0].set_title('Estimada')
axs[0].set_ylabel('Concentracion')
axs[1].plot(x2, pm25, 'tab:red')
axs[1].set_title('Exacta')
axs[1].set_ylabel('Concentracion')
plt.xlabel('Iteraciones')
plt.show()

#################################################################################################
#### CARGAMOS LOS DATOS DE PM
f = lambda x : (x.replace(",","."))
# PM 10
data = pd.read_csv("visual/PM_estimado.csv", sep = ",", encoding = "ISO=8859-1", decimal = ".")
result = data.values.tolist()
pm10 = []

for a in result:
    value = a
    pm10.append(value)


# PM 25
data = pd.read_csv("visual/PM_exacto.csv", sep = ",", encoding = "ISO=8859-1", decimal = ".")
result = data.values.tolist()
pm25 = []

# for a in result:
#     value = a
#     pm25.append(value)

# largo = len(pm10)
# x = np.arange(0, largo)
# plt.plot(x,pm10)
# plt.title("PM 10 - Conservatividad")
# plt.show()

# largo = len(pm25)
# x = np.arange(0, largo)
# plt.plot(x,pm25)
# plt.title("PM 2.5 - Conservatividad")
# plt.show()