import numpy as np
import openmesh as om
import pandas as pd

mesh = om.TriMesh()

# vertices de malla
vh0 = mesh.add_vertex([0, 0, 0])
vh1 = mesh.add_vertex([0, 10, 0])
vh2 = mesh.add_vertex([10, 0, 0])
vh3 = mesh.add_vertex([10, 10, 0])

vh4 = mesh.add_vertex([5, 0, 0])
vh5 = mesh.add_vertex([6, 10, 0])

vh6 = mesh.add_vertex([3, 5, 0])
vh7 = mesh.add_vertex([8, 5, 0])

#Celdas
fh0 = mesh.add_face(vh1, vh6, vh5)
fh1 = mesh.add_face(vh1, vh0, vh6)
fh2 = mesh.add_face(vh0, vh4, vh6)
fh3 = mesh.add_face(vh4, vh7, vh6)
fh4 = mesh.add_face(vh7, vh5, vh6)
fh5 = mesh.add_face(vh7, vh2, vh5)
fh6 = mesh.add_face(vh3, vh2, vh7)
fh7 = mesh.add_face(vh4, vh3, vh7)

# recorremos los lados para contarlos
for eh in mesh.edges():
    nEdges = eh.idx() #devuelve el indice de un lado


numEdges = nEdges + 1 #nuemo total de lados


for fh in mesh.vertices():
    nCell = fh.idx() #devuelve el indice de la celda

    aVetX = [] # guardar las coordenadas X de los vertices
    aVetY = [] # guardar las coordenadas Y de los vertices

    for vh in mesh.fv(fh):
        nVertex = vh.idx()

    print('OK')

    for vh in mesh.fv(fh):
        punto = mesh.point(vh) # devolver un arreglo de long 3
        vetX = punto[0]
        vetY = punto[1]
        aVetX.append(vetX)
        aVetY.append(vetY)
    # aVetX y aVetY deberian tener cada una de las coordenadas de los vertices
    iaux = 0 # contador auxiliar