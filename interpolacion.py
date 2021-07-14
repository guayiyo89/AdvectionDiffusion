import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


aTemp_PLC = []
aTemp_NL = []
aTemp_LE = []

aLluvia_PLC = []
aLluvia_NL = []
aLluvia_LE = []

aVviento_PLC = []
aVviento_NL = []
aVviento_LE = []

aDviento_PLC = []
aDviento_NL = []
aDviento_LE = []

aPresion_PLC = []
aPresion_NL = []
aPresion_LE = []

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

def valEstimado(lim1, lim2, t, t1):
    val = lim1 + (lim2 - lim1) * (t - t1)
    return val

####################################################################################################################
# TEMPERATURA

# Las Encinas
data = pd.read_csv("datos\LE_temp.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
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
    aTemp_LE.append(b)

temp_LE = interpolacion(vacios, aTemp_LE)

# Nielol
data = pd.read_csv("datos\CNL_temp.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
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
    aTemp_NL.append(b)

temp_NL = interpolacion(vacios, aTemp_NL)

# Padre las Casas
data = pd.read_csv("datos\PLC_temp.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
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
    aTemp_PLC.append(b)

temp_PLC = interpolacion(vacios, aTemp_PLC)

####################################################################################################################
# LLUVIA

# Las Encinas
data = pd.read_csv("datos\LE_lluvia.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
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
    aLluvia_LE.append(b)

lluvia_LE = interpolacion(vacios, aLluvia_LE)

# Nielol
data = pd.read_csv("datos\CNL_lluvia.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
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
    aLluvia_NL.append(b)

lluvia_NL = interpolacion(vacios, aLluvia_NL)

# Padre las Casas
data = pd.read_csv("datos\PLC_lluvia.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
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
    aLluvia_PLC.append(b)

lluvia_PLC = interpolacion(vacios, aLluvia_PLC)

####################################################################################################################
# VELC VIENTO

# Las Encinas
data = pd.read_csv("datos\LE_vviento.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
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
    aVviento_LE.append(b)

velViento_LE = interpolacion(vacios, aVviento_LE)

# Nielol
data = pd.read_csv("datos\CNL_vviento.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
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
    aVviento_NL.append(b)

velViento_NL = interpolacion(vacios, aVviento_NL)

# Padre las Casas
data = pd.read_csv("datos\PLC_vviento.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
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
    aVviento_PLC.append(b)

velViento_PLC = interpolacion(vacios, aVviento_PLC)

####################################################################################################################
#  DIR VIENTO

# Las Encinas
data = pd.read_csv("datos\LE_dviento.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
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
    aDviento_LE.append(b)

dirViento_LE = interpolacion(vacios, aDviento_LE)

# Nielol
data = pd.read_csv("datos\CNL_dviento.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
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
    aDviento_NL.append(b)

dirViento_NL = interpolacion(vacios, aDviento_NL)

# Padre las Casas
data = pd.read_csv("datos\PLC_dviento.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
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
    aDviento_PLC.append(b)

dirViento_PLC = interpolacion(vacios, aDviento_PLC)

####################################################################################################################
#  PRESION

# Las Encinas
data = pd.read_csv("datos\LE_presion.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
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
    aPresion_LE.append(b)

presion_LE = interpolacion(vacios, aPresion_LE)

# Nielol
data = pd.read_csv("datos\CNL_presion.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
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
    aPresion_NL.append(b)

presion_NL = interpolacion(vacios, aPresion_NL)

# Padre las Casas
data = pd.read_csv("datos\PLC_presion.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
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
    aPresion_PLC.append(b)

presion_PLC = interpolacion(vacios, aPresion_PLC)

################################################################################################################
# FUNCION

def resultado(cellAs,value):
    if cellAs == 'LE':
        if value == 'temp':
            return temp_LE
        if value == 'rain':
            return lluvia_LE
        if value == 'windVel':
            return velViento_LE
        if value == 'windDir':
            return dirViento_LE
        if value == 'presion':
            return presion_LE
    if cellAs == 'PLC':
        if value == 'temp':
            return temp_PLC
        if value == 'rain':
            return lluvia_PLC
        if value == 'windVel':
            return velViento_PLC
        if value == 'windDir':
            return dirViento_PLC
        if value == 'presion':
            return presion_PLC
    if cellAs == 'NL':
        if value == 'temp':
            return temp_NL
        if value == 'rain':
            return lluvia_NL
        if value == 'windVel':
            return velViento_NL
        if value == 'windDir':
            return dirViento_NL
        if value == 'presion':
            return presion_NL

##########################################################
#GRAFICO DE VIENTO
# largo = len(lluvia)
# x = np.arange(0, largo)
# plt.plot(x,lluvia)
# plt.ylabel('mm3')
# plt.xlabel('Horas')
# plt.title("Precipitaciones")
# plt.show()


ejemplo = resultado('PLC', 'rain')
print(ejemplo)