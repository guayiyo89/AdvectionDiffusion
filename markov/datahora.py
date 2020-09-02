import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Velc viento
data = pd.read_csv(".\markov\pm10_190601_190901_hora.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
result = data.values.tolist()
velv = []

dataPM10 = []
hora10 = []
for a in result:
    if pd.notna(a[2]):
        dataPM10.append(a[2])
        hora10.append(a[1])
    if pd.notna(a[3]):
        dataPM10.append(a[3])
        hora10.append(a[1])
    if pd.notna(a[4]):
        dataPM10.append(a[4])
        hora10.append(a[1])


data = pd.read_csv(".\markov\pm25_190601_190901_hora.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
result1 = data.values.tolist()
velv = []

dataPM25 = []
hora25 = []
for a in result1:
    if pd.notna(a[2]):
        dataPM25.append(a[2])
        hora25.append(a[1])
    if pd.notna(a[3]):
        dataPM25.append(a[3])
        hora25.append(a[1])
    if pd.notna(a[4]):
        dataPM25.append(a[4])
        hora25.append(a[1])
  

def calcular(data,hora):
    time = np.linspace(0,2300,24)

    promedios = []

    for t in time:
        values = []
        for i in range(len(hora)):
            if t == hora[i]:
                values.append(data[i])
        prom = sum(values) / len(values)
        promedios.append(prom)
    
    return promedios

promPM10 = calcular(dataPM10,hora10)
promPM25 = calcular(dataPM25,hora25)

horas = np.linspace(0,23,24)

fig, axs = plt.subplots(2, sharex = True)
fig.suptitle('Promedios Concentraciones')
axs[0].plot(horas, promPM10)
axs[0].set_title('PM 10')
axs[0].set_ylabel('Concentracion (ug/m3)')
axs[1].plot(horas, promPM25, 'tab:red')
axs[1].set_title('PM 2.5')
axs[1].set_ylabel('Concentracion (ug/m3)')
plt.xlabel('Hora')
plt.show()

print('PM2.5: ', promPM25)
print('PM10: ', promPM10)