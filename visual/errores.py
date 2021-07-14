import numpy as np
from numpy.core.fromnumeric import mean
import pandas as pd
import matplotlib.pyplot as plt



#### CARGAMOS LOS DATOS DE PM
f = lambda x : (x.replace(",","."))
# PM 10
data = pd.read_csv(r"visual/PM_estimadoExb.csv", sep = ",", encoding = "ISO=8859-1", decimal = ".")
result = data.values.tolist()
estimado = []

data = pd.read_csv(r"visual/PM_estimadoImp.csv", sep = ",", encoding = "ISO=8859-1", decimal = ".")
result1 = data.values.tolist()
estimado5 = []

cornerN = [17,5]
cornerS = [27,10]
borderS = [43,45,48,71]
borderN = [51,65,38,78]
borderO = [21,24,69]
borderE = [57,34,77]

interiorO = [966,491,819,135]
interiorE = [656,689,668,1502]
centers = [1327, 812, 1100]

cells = cornerN + cornerS

# PM 25
data1 = pd.read_csv(r"visual/PM_exactoExb.csv", sep = ",", encoding = "ISO=8859-1", decimal = ".")
result2 = data1.values.tolist()
exacto = []

data1 = pd.read_csv(r"visual/PM_exactoImp.csv", sep = ",", encoding = "ISO=8859-1", decimal = ".")
result3 = data1.values.tolist()
exacto5 = []

aCell_error = [] #guardo errores porcentuales promedios por cada celda
aCell_error2 = [] #guardo errores porcentuales promedios por cada celda

for cell in cells:
    estimado = []
    exacto = []
    estimado5 = []
    exacto5 = []
    aErr_ex = []
    aErr_im = []

    for a in result:
        value = a[cell]
        estimado.append(value)

    for a in result2:
        value = a[cell]
        exacto.append(value)

    ## obtengo valores de t (Parar t = 2.5)
    aTmp = np.linspace(0, 3.0, len(estimado))
    for i in range(len(exacto)):
        #error relativo porcentual
        err = (abs(exacto[i] - estimado[i]) / exacto[i])*100
        aErr_ex.append(err)

    aCell_error.append(mean(aErr_ex))
    
    for a in result1:
        value = a[cell]
        estimado5.append(value)

    for a in result3:
        value = a[cell]
        exacto5.append(value)

    ## obtengo valores de t (Parar t = 2.5)
    aTmp = np.linspace(0, 3.0, len(estimado))

    for i in range(len(exacto5)):
        #error relativo porcentual
        err = (abs(exacto5[i] - estimado5[i]) / exacto5[i])*100
        aErr_im.append(err)

    aCell_error2.append(mean(aErr_im))

#################################################################################################
print('Error relativo porcentual')
print('explicito: ',aCell_error)
print('implicito: ',aCell_error2)