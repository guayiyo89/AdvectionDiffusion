import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# data PM2.5
data = pd.read_csv("PM_Estimado.csv")
result = data.values.tolist()
datanew = np.asarray(result)

dataPM25 = []

for a in result:
    c = -1 # No hay lluvia negativa
    if pd.notna(a[12]):
        c = a[12]

    b = c
    dataPM25.append(b)

m = np.size(dataPM25)

dataPM25 = np.array(dataPM25)
# dataPM25 = dataPM25 * pow(10,9)

# data PM10
data = pd.read_csv("PM_exacto.csv")
result = data.values.tolist()
datanew = np.asarray(result)

dataPM10 = []

for a in result:
    c = -1 # No hay lluvia negativa
    if pd.notna(a[12]):
        c = a[12]

    b = c
    dataPM10.append(b)


# data Temperatura
data = pd.read_csv("temperatura_test.csv")
result = data.values.tolist()
datanew = np.asarray(result)

dataTmp = []

for a in result:
    valor = a[0]
    dataTmp.append(valor)


# data Viento
data = pd.read_csv("viento_test.csv")
result = data.values.tolist()
datanew = np.asarray(result)

dataViento = []

for a in result:
    valor = a[0]
    dataViento.append(valor)


dataPM10 = np.array(dataPM10)

print(len(dataTmp))
print(len(dataViento))

x = np.arange(1,m+1)

fig, axs = plt.subplots(2)
fig.suptitle('Vertically stacked subplots')
axs[0].plot(x,dataPM25)
axs[0].set_ylabel('PM 2.5 (ug/m^3)')
plt.ylabel('PM 2.5 (ug/m^3)')
axs[1].plot(x,dataPM10, 'tab:red')
axs[1].set_ylabel('PM 10 (ug/m^3)')
axs[1].set_xlabel('Tiempo (horas)')

plt.show()


fig, axs = plt.subplots(2)
fig.suptitle('Vertically stacked subplots')
axs[0].plot(x,dataViento, 'tab:green')
axs[0].set_ylabel('Viento (m/h)')
axs[1].plot(x,dataTmp)
axs[1].set_ylabel('Temperatura (C)')
axs[1].set_xlabel('Tiempo (horas)')

plt.show()