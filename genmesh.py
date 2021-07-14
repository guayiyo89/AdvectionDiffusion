import sys

from numpy.core.shape_base import vstack
import openmesh as om
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile


def save_mesh(P, T, filename):
    mesh = om.TriMesh()
    
    filename2 = filename[0:-4]+'_2.off'
    
    print('WARNING: All prints are printed into a file.')
    sys.stdout = open(filename2, 'w')
    print('OFF')
    print(len(P), len(T.simplices), 0)
    
    k=0
    vh=[]
    # de un arreglo de puntos genera los vertices
    for p in P:
        k=k+1; 
        print('%.2f' % p[0] ,'%.2f' % p[1], 0)    
        vh.append(mesh.add_vertex([p[0], p[1], 0]))
    
    k=0
    # generamos las caras
    for t in T.simplices:
        k=k+1
        print(3,' ',t[0], t[1], t[2])
        mesh.add_face(vh[t[0]],vh[t[1]],vh[t[2]])

    om.write_mesh(filename, mesh)

X = np.arange(-5,5.5,0.5)
Y = np.arange(-5,5.5,0.5)

Xbar = np.arange(-4.75,5,0.5)
Ybar = np.arange(-4.75,5,0.5)
print(X,Y)
print(Xbar, Ybar)

# X = np.arange(0,5,0.1)
# Y = np.arange(0,1.05,0.02)

P = np.array([[x,y] for y in Y for x in X])
Pbar = np.array([[x,y] for y in Ybar for x in Xbar])
PP = vstack([P, Pbar])
T = Delaunay(PP)
plt.triplot(PP[:,0],PP[:,1],T.simplices)
#plt.plot(P[:,0],P[:,1],'go')

save_mesh(PP, T, 'pruebacell.off')

mesh = om.read_trimesh('prueba1000cell.off')

#CALCULAR LONGITUD DE LOS LADOS DE LAS CELDAS
Xbar =[]
Ybar = []

#Recorrer todas las interfaces
for eh in mesh.edges():
    totEd = eh.idx()

totEd = totEd + 1

Xe  = np.zeros(totEd)             # arreglos que guardan los centros de cada lado
Ye  = np.zeros(totEd)

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
    Xbar.append(barX)
    Ybar.append(barY)

    iaux = 0 

    for eh in mesh.fe(fh):
        nEdge = eh.idx()
        nCentLadoX = 0.5 * (aVetX[iaux - 1] + aVetX[iaux])
        nCentLadoY = 0.5 * (aVetY[iaux - 1] + aVetY[iaux])
        Xe[nEdge] = nCentLadoX
        Ye[nEdge] = nCentLadoY
        iaux = iaux + 1


##### DIBUJAR

for fh in mesh.faces():   
    nCell = fh.idx()
    plt.plot(Xbar[nCell], Ybar[nCell], marker = 'o', color ='thistle', markersize=15)
    #plt.text(Xbar[nCell], Ybar[nCell], nCell, fontsize=11, horizontalalignment='center', verticalalignment='center')

for vh in mesh.vertices():
    p = mesh.point(vh)
    
    # print(eh.idx(), mesh.is_boundary(vh), point_is_front[vh.idx()]) 

    # if point_is_front[vh.idx()]==1:
    if mesh.is_boundary(vh):    
            color = 'orchid'
            markersize=8
    else:
            color = 'lavender'
            markersize=8
            
    #plt.plot(p[0], p[1], 'o', markersize=markersize, color='green')
    #plt.text(p[0], p[1], vh.idx(), fontsize=12, horizontalalignment='center', verticalalignment='center')

for eh in mesh.edges():
    nEdge = eh.idx()
    # if edge_is_front[nEdge]==1:
    if mesh.is_boundary(eh):
        color = 'orange'
    else:
        color = 'yellow'
        
    #plt.plot(Xe[nEdge], Ye[nEdge], marker = 'o', color = color, markersize=10)
    #plt.text(Xe[nEdge], Ye[nEdge], nEdge, fontsize=10, horizontalalignment='center', verticalalignment='center')

plt.show()