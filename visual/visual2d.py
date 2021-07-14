import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import openmesh as om
from scipy.spatial import Delaunay

#### CARGAMOS LOS DATOS DE PM
f = lambda x : (x.replace(",","."))
# PM 10
data = pd.read_csv(r"visual/PM_estimadoExpT.csv", sep = ",", encoding = "ISO=8859-1", decimal = ".")
result = data.values.tolist()
estimado = []

data = pd.read_csv(r"visual/PM_exactoExpT.csv", sep = ",", encoding = "ISO=8859-1", decimal = ".")
result1 = data.values.tolist()
exacto = []


########################################
### Cargo la Malla
mesh = om.read_trimesh('prueba1000cell.off')

#baricentros
Xbars = []
Ybars = []

for fh in mesh.faces():
    nCell = fh.idx()

    #guardo las coordenadas de los vertices
    aVetX = []
    aVetY = []

    for vh in mesh.fv(fh):
        punto = mesh.point(vh)
        vetX = punto[0]
        vetY = punto[1]
        aVetX.append(vetX)
        aVetY.append(vetY)

    # Calcular el baricentro de la celda
    barX = sum(aVetX)/3
    barY = sum(aVetY)/3
    Xbars.append(barX)
    Ybars.append(barY)

    iaux = 0 

X = np.arange(-5,5.5,0.5)
Y = np.arange(-5,5.5,0.5)

Xbar = np.arange(-4.75,5,0.5)
Ybar = np.arange(-4.75,5,0.5)

# X = np.arange(0,5,0.1)
# Y = np.arange(0,1.05,0.02)


P = np.array([[x,y] for y in Y for x in X])
Pbar = np.array([[x,y] for y in Ybar for x in Xbar])
PP = np.vstack([P, Pbar])
T = Delaunay(PP)

print(len(result))

aTimes = [0, 200, 400, 600, 800, 1000, 1250, 1500, 1750, 2000, 2174]
#aTimes = [2000, 2174]

tiempos = np.linspace(0,3,len(result))

for t in aTimes:
    uPM_est = result[t]
    uPM_ex = result1[t]
    tVal = tiempos[t]

    fig, axs = plt.subplots(2, sharex = True)
    fig.suptitle(f'Explicito CFL=0.75 (Diff 0.001) - t = {tVal}')


    for fh in mesh.faces():   
        nCell = fh.idx()
        val = uPM_est[nCell]
        if val <= 1.025:
            col = '#137e48'
        elif val > 1.025 and val <= 1.05:
            col = '#329040'
        elif val > 1.05 and val <= 1.1:
            col = '#53a239'
        elif val > 1.1 and val <= 1.2:
            col = '#87bf2d'
        elif val > 1.2 and val <= 1.3:
            col = '#bcd821'
        elif val > 1.3 and val <= 1.4:
            col = '#e2cb18'
        elif val > 1.4 and val <= 1.5:
            col = '#e0e018'
        elif val > 1.5 and val <= 1.6:
            col = '#e2ce18'
        elif val > 1.6 and val <= 1.7:
            col = '#e2b718'
        elif val > 1.7 and val <= 1.8:
            col = '#e2a318'
        elif val > 1.8:
            col = '#e28a18'

        axs[0].plot(Xbars[nCell], Ybars[nCell], marker = 'o', color = col, markersize=7.6)
        val = uPM_ex[nCell]
        if val <= 1.025:
            col = '#137e48'
        elif val > 1.025 and val <= 1.05:
            col = '#329040'
        elif val > 1.05 and val <= 1.1:
            col = '#53a239'
        elif val > 1.1 and val <= 1.2:
            col = '#87bf2d'
        elif val > 1.2 and val <= 1.3:
            col = '#bcd821'
        elif val > 1.3 and val <= 1.4:
            col = '#e2cb18'
        elif val > 1.4 and val <= 1.5:
            col = '#e0e018'
        elif val > 1.5 and val <= 1.6:
            col = '#e2ce18'
        elif val > 1.6 and val <= 1.7:
            col = '#e2b718'
        elif val > 1.7 and val <= 1.8:
            col = '#e2a318'
        else:
            col = '#e28a18'

        axs[1].plot(Xbars[nCell], Ybars[nCell], marker = 'o', color = col, markersize=7.6)
        #plt.text(Xbars[nCell], Ybars[nCell], nCell, fontsize=11, horizontalalignment='center', verticalalignment='center')


    axs[0].triplot(PP[:,0],PP[:,1],T.simplices, color= 'white')
    axs[0].set_title('Estimado')
    axs[1].triplot(PP[:,0],PP[:,1],T.simplices, color= 'white')
    axs[1].set_title('Exacto')

    plt.show()