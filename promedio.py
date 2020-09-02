import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#### CARGAMOS LOS DATOS DE PM
f = lambda x : (x.replace(",","."))
# PM 2.5
data = pd.read_csv("pmdatos\dnielolDia_25.csv", sep = ",", encoding = "ISO=8859-1")
result = data.values.tolist()
pm25 = []

for a in result:
    c = -1 # default
    if pd.notna(a[1]):
        c = a[1]

    b = c
    pm25.append(b)

# PM 10
data = pd.read_csv("pmdatos\dnielolDia_10.csv", sep = ",", encoding = "ISO=8859-1")
result = data.values.tolist()
pm10 = []

for a in result:
    c = -1 # default
    if pd.notna(a[1]):
        c = a[1]

    b = c
    pm10.append(b)


PM25_f = []
PM10_f = []

for value in pm10:
    if value != -1:
        PM10_f.append(value)

for value in pm25:
    if value != -1:
        PM25_f.append(value)

Prom25 = np.mean(PM25_f)
Prom10 = np.mean(PM10_f)

#PM25_f = PM25_f.sort()

Med25 = np.median(PM25_f)
Med10 = np.median(PM10_f)
print(Prom10)
print(Prom25)
print("---------------")
print(Med10)
print(Med25)

