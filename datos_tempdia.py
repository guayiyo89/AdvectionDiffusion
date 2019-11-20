import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#### CARGAMOS LOS DATOS DE PM
f = lambda x : (x.replace(",","."))

# data
data = pd.read_csv("datos\datos_lluvia_encinas.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
result = data.values.tolist()
datanew = np.asarray(result)

data = []

for a in result:
    c = -1 # No hay lluvia negativa
    if pd.notna(a[2]):
        c = a[2]

    b = c
    data.append(b)

m = np.size(data)

temp_dia = []
auxiliar = []
i = 0
for valor in data:
    if i == 24:
        prom = np.mean(auxiliar)
        temp_dia.append(prom)
        auxiliar = []
        i = 0

    auxiliar.append(valor)

    i = i + 1

df2 = pd.DataFrame(temp_dia, columns = ['Registros'])

export = df2.to_csv(r'lluvia_Encina_dia.csv', index = None, header=True)

print(df2)