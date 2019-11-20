import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Direction viento
data = pd.read_csv("datos\datos_temp_encinas.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
result = data.values.tolist()
temperatura = []

for a in result:
    if pd.notna(a[2]):
        temperatura.append(a[2])

m = np.size(temperatura)
print(m)

temp = np.asarray(temperatura)
print('MAX y MIN')
print(max(temp)); print(min(temp))
print('PROM y MEDIAN')
print(np.mean(temp)); print(np.median(temp))
#PROM = 11.39 y MEDIAN = 11.1

muybaja = []; baja = []; mediabaja = []; media = []; mediaalta = []; alta = []; muyalta = []; extrema = []

for t in temp:
    if t < 0:
        muybaja.append(t)
    if t >=0 and t < 5:
        baja.append(t)
    if t >=5 and t < 10:
        mediabaja.append(t)
    if t >=10 and t < 15:
        media.append(t)
    if t >=15 and t < 20:
        mediaalta.append(t)
    if t >=20 and t < 25:
        alta.append(t)
    if t >=25 and t <= 30:
        muyalta.append(t)
    if t >30:
        extrema.append(t)

print(len(muybaja)); print(len(baja)); print(len(mediabaja)); print(len(media))
print(len(mediaalta)); print(len(alta)); print(len(muyalta)); print(len(extrema))