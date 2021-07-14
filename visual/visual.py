import numpy as np
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

cornerN = [17,18,5,6]
cornerS = [27,10]
borderS = [43,45,48,71]
borderN = [51,65,38,78]
borderO = [21,24,69]
borderE = [57,34,77]

interiorO = [966,491,819,135,222]
interiorE = [656,689,668,1502,1450]
centers = [1327, 812, 1220,1100]

cells = centers

# PM 25
data1 = pd.read_csv(r"visual/PM_exactoExb.csv", sep = ",", encoding = "ISO=8859-1", decimal = ".")
result2 = data1.values.tolist()
exacto = []

data1 = pd.read_csv(r"visual/PM_exactoImp.csv", sep = ",", encoding = "ISO=8859-1", decimal = ".")
result3 = data1.values.tolist()
exacto5 = []

for cell in cells:
    estimado = []
    exacto = []
    estimado5 = []
    exacto5 = []

    for a in result:
        value = a[cell]
        estimado.append(value)

    for a in result2:
        value = a[cell]
        exacto.append(value)

    for a in result1:
        value = a[cell]
        estimado5.append(value)

    for a in result3:
        value = a[cell]
        exacto5.append(value)

    largo = len(estimado)
    x = np.arange(0, largo)
    x = np.linspace(0,3,largo)

    largo2 = len(estimado5)
    x2 = np.arange(0, largo2)
    x2 = np.linspace(0,3,largo2)


    fig, axs = plt.subplots(2, sharex = True)
    fig.suptitle(f'Solucion Manufacturada (Diff 0.001) - Celda {cell}')
    axs[0].plot(x, estimado, label='Estimado')
    axs[0].plot(x, exacto, 'tab:red', label='Exacto')
    axs[0].set_title('Explicita CFL 0.75')
    axs[0].yaxis.set_label_position("right")
    axs[0].set_ylabel('Concentracion')
    axs[1].plot(x2, estimado5, label='Estimado')
    axs[1].plot(x2, exacto5, 'tab:red', label='Exacto')
    axs[1].set_title('Implicita')
    axs[1].yaxis.set_label_position("right")
    axs[1].set_ylabel('Concentracion')
    plt.xlabel('tiempo')
    plt.show()

#################################################################################################