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


# Velc viento
data = pd.read_csv("datos\datos_lluvia_encinas.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
result = data.values.tolist()
lluvia = []

for a in result:
    if pd.notna(a[2]):
        lluvia.append(a[2])

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

grupo1 = []; grupo2 = []; grupo3 = []; grupo4 = []; grupo5 = []; grupo6 = []; grupo7 = []; grupo8 = []
grupo9 = []; grupo10 = []; grupo11 = []; grupo12 = []; grupo13 = []; grupo14 = []; grupo15 = []; grupo16 = []
grupo17 = []; grupo18 = []; grupo19 = []; grupo20 = []; grupo21 = []; grupo22 = []; grupo23 = []; grupo24 = []

for valor in dirv:

    if valor >= 0 and valor < 45:
        grupo1.append(valor)
    if valor >= 45 and valor < 90:
        grupo2.append(valor)
    if valor >= 90 and valor < 135:
        grupo3.append(valor)
    if valor >= 135 and valor < 180:
        grupo4.append(valor)
    if valor >= 180 and valor < 225:
        grupo5.append(valor)
    if valor >= 225 and valor < 270:
        grupo6.append(valor)
    if valor >= 270 and valor < 315:
        grupo7.append(valor)
    if valor >= 315 and valor < 360:
        grupo8.append(valor)


print(len(velv))

largos1 = [len(grupo1), len(grupo2), len(grupo3), len(grupo4), len(grupo5), len(grupo6), len(grupo7), len(grupo8)]

print('Largo lluvia: ',len(lluvia))
print('Max lluvia: ',max(lluvia))
print('Min lluvia: ',min(lluvia))
print('Prom lluvia: ',np.mean(lluvia))
print('Mediana lluvia: ',np.median(lluvia))


# dirViento = np.concatenate((grupo11, grupo12, grupo13, grupo14))

sec1 = []; sec2 = []; sec3 = []; sec4 = []; sec5 = []; sec6 = []

for pp in lluvia:
    if pp == 0:
        sec1.append(pp)
    if pp > 0 and pp <= 0.1:
        sec2.append(pp)
    if pp > 0.1 and pp <= 0.25:
        sec3.append(pp)
    if pp > 0.25 and pp <= 1.0:
        sec4.append(pp)
    if pp > 1.0 and pp <= 4.0:
        sec5.append(pp)
    if pp > 4.0:
        sec6.append(pp)

pps = [len(sec1), len(sec2), len(sec3), len(sec4), len(sec5), len(sec6)]

print(pps)

print