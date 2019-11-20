import numpy as np
import pandas as pd

# Direction viento
data = pd.read_csv("datos\datos_dirviento.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
result = data.values.tolist()
dirv = []

for a in result:
    if pd.notna(a[2]):
        dirv.append(a[2])


# Velc viento
data = pd.read_csv("datos\datos_velviento_encinas.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
result = data.values.tolist()
velv = []

for a in result:
    if pd.notna(a[2]):
        velv.append(a[2])


PromVel = np.mean(velv)
PromDir = np.mean(dirv)

MedDir = np.median(dirv)
MedVel = np.median(velv)

print(PromDir)
print(PromVel)
print('---------------')
print(MedDir)
print(MedVel)

################ CLASIFICACION DE LOS VIENTOS ################

grupo1 = []; grupo2 = []; grupo3 = []; grupo4 = []; grupo5 = []; grupo6 = []; grupo7 = []
grupo8 = []
for valor in velv:

    valor = valor * 3.6
    if valor >= 0 and valor < 3:
        grupo1.append(valor)
    if valor >= 3 and valor < 6:
        grupo2.append(valor)
    if valor >= 6 and valor < 9:
        grupo3.append(valor)
    if valor >= 9 and valor < 12:
        grupo4.append(valor)
    if valor >= 12 and valor < 15:
        grupo5.append(valor)
    if valor >= 15 and valor < 20:
        grupo6.append(valor)
    if valor >= 20:
        grupo7.append(valor)

print(len(velv))

largos = [len(grupo1), len(grupo2), len(grupo3), len(grupo4), len(grupo5), len(grupo6), len(grupo7)]

print(largos)

# dirViento = np.concatenate((grupo11, grupo12, grupo13, grupo14))

print(np.mean(velv))
print(np.median(velv))