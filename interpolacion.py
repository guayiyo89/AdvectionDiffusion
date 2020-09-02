import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


aTemp = []
aLluvia = []
aVviento = []
aDviento = []

def interpolacion(emptyvalues,aValues):
    for a in emptyvalues:
        if a == 0 or a == len(aValues):
            continue
        else:
            # limite inferior de la interpolacion
            x0 =  a - 1
            lim_min = aValues[x0]
            #print(lim_min)
            # limite mayor
            x1 = a + 1
            while aValues[x1] == -1:
                x1 = x1 + 1
            lim_max = aValues[x1]
            #print(lim_max)
        
        aValues[a] = lim_min + ((lim_max - lim_min)/(x1 - x0)) * (a - x0)
    
    return aValues

####################################################################################################################
# TEMPERATURA
data = pd.read_csv("ambientales_var\Temp_inv_19.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
result = data.values.tolist()
datanew = np.asarray(result)

vacios = [] #guardo los indices de celdas nulas
i = 0 #index
for a in result:
    c = -276 # Temperatura imposible de llegar
    if pd.notna(a[2]):
        c = a[2]
    else:
        vacios.append(i)

    b = c
    aTemp.append(b)

temperatura = interpolacion(vacios, aTemp)

####################################################################################################################
# LLUVIA
data = pd.read_csv("ambientales_var\lluvia_inv_19.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
result = data.values.tolist()
datanew = np.asarray(result)

vacios = [] #guardo los indices de celdas nulas
i = 0 #index
for a in result:
    c = -1 # No hay lluvia negativa
    if pd.notna(a[2]):
        c = a[2]
    else:
        vacios.append(i)
    b = c
    aLluvia.append(b)

lluvia = interpolacion(vacios, aLluvia)

####################################################################################################################
# VELC VIENTO
data = pd.read_csv("ambientales_var\Viento_inv_19.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
result = data.values.tolist()
datanew = np.asarray(result)

vacios = [] #guardo los indices de celdas nulas
i = 0 #index
for a in result:
    c = -1 # No hay velc negativa
    if pd.notna(a[2]):
        c = a[2]
    else:
        vacios.append(i)
    b = c
    i = i + 1
    aVviento.append(b)

velViento = interpolacion(vacios, aVviento)

####################################################################################################################
#  DIR VIENTO
data = pd.read_csv("ambientales_var\Vientod_inv_19.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
result = data.values.tolist()
datanew = np.asarray(result)

vacios = [] #guardo los indices de celdas nulas
i = 0 #index

for a in result:
    c = -1 # No hay direccion negativa
    if pd.notna(a[2]):
        c = a[2]
    else:
        vacios.append(i)
    b = c
    i = i + 1
    aDviento.append(b)

dirViento = interpolacion(vacios, aDviento)

def resultado(value):
    if value == 'temp':
        return temperatura
    if value == 'rain':
        return lluvia
    if value == 'windVel':
        return velViento
    if value == 'windDir':
        return dirViento

##########################################################
#GRAFICO DE VIENTO
# largo = len(lluvia)
# x = np.arange(0, largo)
# plt.plot(x,lluvia)
# plt.ylabel('mm3')
# plt.xlabel('Horas')
# plt.title("Precipitaciones")
# plt.show()