import numpy as np
import openmesh as om
import pandas as pd
import exactas as fnc

#### Defino el tipo de Malla
mesh = om.TriMesh() # TriMesh es una funcion de OpenMesh (om)


# vertices de malla
vh0 = mesh.add_vertex([-3.8602, 5.0, 0])
vh1 = mesh.add_vertex([-3.8602, -3.1369, 0])
vh2 = mesh.add_vertex([3.8602, -3.1369, 0])
vh3 = mesh.add_vertex([3.8602, 5.0, 0])

vh4 = mesh.add_vertex([-2.1695, 5.0, 0])
vh5 = mesh.add_vertex([-0.0958, 5.0, 0])
vh6 = mesh.add_vertex([2.1839, 5.0, 0])
vh7 = mesh.add_vertex([-2.1695, -3.1369, 0])
vh8 = mesh.add_vertex([-0.0958, -3.1369, 0])
vh9 = mesh.add_vertex([2.1839, -3.1369, 0])

vh10 = mesh.add_vertex([-3.8602, 3.1290, 0])
vh11 = mesh.add_vertex([-3.8602, 0.8997, 0])
vh12 = mesh.add_vertex([-3.8602, -0.9236, 0])
vh13 = mesh.add_vertex([3.8602, 3.1290, 0])
vh14 = mesh.add_vertex([3.8602, 0.8997, 0])
vh15 = mesh.add_vertex([3.8602, -0.9236, 0])

vh16 = mesh.add_vertex([-2.4234, 3.551, 0])
vh17 = mesh.add_vertex([-0.6561, 3.4395, 0])
vh18 = mesh.add_vertex([1.0105, 3.6385, 0])
vh19 = mesh.add_vertex([2.3994, 3.3519, 0])

vh20 = mesh.add_vertex([-1.887, 1.7914, 0])
vh21 = mesh.add_vertex([0.0766, 1.3933, 0])
vh22 = mesh.add_vertex([1.3745, 1.9268, 0])
vh23 = mesh.add_vertex([2.4425, 1.1306, 0])

vh24 = mesh.add_vertex([-2.7203, -0.6768, 0])
vh25 = mesh.add_vertex([-1.1398, 0.1513, 0])
vh26 = mesh.add_vertex([-0.7902, -1.2261, 0])
vh27 = mesh.add_vertex([0.9435, -0.6051, 0])
vh28 = mesh.add_vertex([2.3994, -0.6927, 0])

vh29 = mesh.add_vertex([-5.0, 5.0, 0])
vh30 = mesh.add_vertex([-5.0, -5.0, 0])
vh31 = mesh.add_vertex([5.0, -5.0, 0])
vh32 = mesh.add_vertex([5.0, 5.0, 0])

vh33 = mesh.add_vertex([-5.0, 2.3328, 0])
vh34 = mesh.add_vertex([-5.0, 0.0557, 0])
vh35 = mesh.add_vertex([-5.0, -1.9268, 0])
vh36 = mesh.add_vertex([5.0, 2.3328, 0])
vh37 = mesh.add_vertex([5.0, 0.0557, 0])
vh38 = mesh.add_vertex([5.0, -1.9268, 0])

vh39 = mesh.add_vertex([-3.8602, -5.0, 0])
vh40 = mesh.add_vertex([-2.1695, -5.0, 0])
vh41 = mesh.add_vertex([-0.0958, -5.0, 0])
vh42 = mesh.add_vertex([2.1839, -5.0, 0])
vh43 = mesh.add_vertex([3.8602, -5.0, 0])

#Celdas
fh0 = mesh.add_face(vh0, vh10, vh16)
fh1 = mesh.add_face(vh0, vh16, vh4)
fh2 = mesh.add_face(vh4, vh16, vh17)
fh3 = mesh.add_face(vh4, vh17, vh5)
fh4 = mesh.add_face(vh5, vh17, vh18)
fh5 = mesh.add_face(vh5, vh18, vh6)
fh6 = mesh.add_face(vh6, vh18, vh19)
fh7 = mesh.add_face(vh6, vh19, vh3)
fh8 = mesh.add_face(vh3, vh19, vh13)

fh9 = mesh.add_face(vh16, vh10, vh20)
fh10 = mesh.add_face(vh17, vh16, vh20)
fh11 = mesh.add_face(vh17, vh20, vh21)
fh12 = mesh.add_face(vh17, vh21, vh18)
fh13 = mesh.add_face(vh18, vh21, vh22)
fh14 = mesh.add_face(vh18, vh22, vh19)
fh15 = mesh.add_face(vh23, vh19, vh22)
fh16 = mesh.add_face(vh13, vh19, vh23)
fh17 = mesh.add_face(vh14, vh13, vh23)

fh18 = mesh.add_face(vh20, vh10, vh11)
fh19 = mesh.add_face(vh20, vh11, vh24)
fh20 = mesh.add_face(vh20, vh24, vh25)
fh21 = mesh.add_face(vh20, vh25, vh21)
fh22 = mesh.add_face(vh21, vh25, vh27)
fh23 = mesh.add_face(vh27, vh22, vh21)
fh24 = mesh.add_face(vh22, vh27, vh23)
fh25 = mesh.add_face(vh23, vh27, vh28)
fh26 = mesh.add_face(vh23, vh28, vh14)
fh27 = mesh.add_face(vh14, vh28, vh15)

fh28 = mesh.add_face(vh11, vh12, vh24)
fh29 = mesh.add_face(vh26, vh25, vh24)
fh30 = mesh.add_face(vh27, vh25, vh26)

fh31 = mesh.add_face(vh12, vh1, vh24)
fh32 = mesh.add_face(vh24, vh1, vh7)
fh33 = mesh.add_face(vh26, vh24, vh7)
fh34 = mesh.add_face(vh26, vh7, vh8)
fh35 = mesh.add_face(vh26, vh8, vh27)
fh36 = mesh.add_face(vh9, vh27, vh8)
fh37 = mesh.add_face(vh27, vh9, vh28)
fh38 = mesh.add_face(vh28, vh9, vh2)
fh39 = mesh.add_face(vh15, vh28, vh2)

fh40 = mesh.add_face(vh0, vh29, vh10)
fh41 = mesh.add_face(vh29, vh33, vh10)
fh42 = mesh.add_face(vh33, vh11, vh10)
fh43 = mesh.add_face(vh33, vh34, vh11)
fh44 = mesh.add_face(vh11, vh34, vh12)
fh45 = mesh.add_face(vh34, vh35, vh12)
fh46 = mesh.add_face(vh12, vh35, vh1)
fh47 = mesh.add_face(vh35, vh30, vh1)
fh48 = mesh.add_face(vh30, vh39, vh1)

fh49 = mesh.add_face(vh39, vh40, vh1)
fh50 = mesh.add_face(vh40, vh7, vh1)
fh51 = mesh.add_face(vh40, vh8, vh7)
fh52 = mesh.add_face(vh40, vh41, vh8)
fh53 = mesh.add_face(vh41, vh42, vh8)
fh54 = mesh.add_face(vh42, vh9, vh8)
fh55 = mesh.add_face(vh42, vh2, vh9)
fh56 = mesh.add_face(vh42, vh43, vh2)

fh57 = mesh.add_face(vh43, vh31, vh2)
fh58 = mesh.add_face(vh31, vh38, vh2)
fh59 = mesh.add_face(vh38, vh15, vh2)
fh60 = mesh.add_face(vh38, vh37, vh15)
fh61 = mesh.add_face(vh15, vh37, vh14)
fh62 = mesh.add_face(vh37, vh36, vh14)
fh63 = mesh.add_face(vh36, vh13, vh14)
fh64 = mesh.add_face(vh36, vh32, vh13)
fh65 = mesh.add_face(vh32, vh3, vh13)

########################################################################################
#Contar vertices, faces, edges

#Recorro todo los nodo de la malla
for vh in mesh.vertices():
    totVx = vh.idx()

totVx = totVx + 1

#Recorrer todas las interfaces
for eh in mesh.edges():
    totEd = eh.idx()

totEd = totEd + 1

for fh in mesh.faces():
    totCell = fh.idx()

totCell = totCell + 1

########################################################################################
# Variables Geometricas

aLongLados = np.zeros(totEd)             # arreglo q guardaba la longitud de los lados

aCenint_x  = np.zeros(totEd)             # arreglos que guardan los centros de cada lado
aCenint_y  = np.zeros(totEd)

aBarX = []                               # arreglo guarda los baricentros de la celda
aBarY = []

aArea = []                               # arreglo q guarda las areas

aNormal_xf = []                          # arreglos q guardasn las normales
aNormal_yf = []

aDist_cell_edge = []                     # Distancia centro celda a centro lado

aSum_dist_node_cell = []                 # sumatoria de los inversos de las distancias nodo - centro (vertices)


aDist_node_cell = np.zeros((totVx,totCell))             # distancia nodo - centro de celda

link_cell_to_edge = np.zeros((totCell,3))               # relacion celda a sus bordes
link_edge_to_cell = []

link_cell_to_cell = np.zeros((totCell,3))               # relacion celda a celda vecinas

link_cell_to_vertex = np.zeros((totCell,3))             # relacion celda vertices

ladoPar = np.zeros(totCell) - 1                         # paridad lado frontera TEST
celdaPar = np.zeros(totCell) - 1                        # paridad celda frontera TEST

############################################################################################################
############################################################################################################
# INFO GEOMETRICA
############################################################################################################

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
    aBarX.append(barX)
    aBarY.append(barY)

    # Area de la celda
    nArea = aVetX[0] * (aVetY[2] - aVetY[1]) + aVetX[1] * (aVetY[0] - aVetY[2]) + aVetX[2] * (aVetY[1] - aVetY[0])
    nArea = 0.5 * np.absolute(nArea)
    aArea.append(nArea)

    # Variables
    aux_dist_cCell_cEdge = []               # auxiliar para el calculo distancia centro celda - centro lado
    iaux = 0                                # contador auxiliar
    auxNormal_x = []                        # auxiliar q guarda las normales (largo 3)
    auxNormal_y = []

    for eh in mesh.fe(fh):
        nEdge = eh.idx()

        # Calcular longitud
        distX = aVetX[iaux -1] - aVetX[iaux]
        distY = aVetY[iaux -1] - aVetY[iaux]
        nLong = np.sqrt(pow(distX,2) + pow(distY,2))
        aLongLados[nEdge] = nLong                           # guardamos en el arreglo

        # calculado el punto medio de cada lado
        nCentLadoX = 0.5 * (aVetX[iaux - 1] + aVetX[iaux])
        nCentLadoY = 0.5 * (aVetY[iaux - 1] + aVetY[iaux])
        #lo guardamos en el arreglo
        aCenint_x[nEdge] = nCentLadoX
        aCenint_y[nEdge] = nCentLadoY

        # Calcular la distancia baricentro a centro de segmento
        dif_x = aBarX[nCell] - nCentLadoX
        dif_y = aBarY[nCell] - nCentLadoY
        distanceCL = np.sqrt(pow(dif_x,2) + pow(dif_y,2))
        aux_dist_cCell_cEdge.append(distanceCL)             # se guarda

        # tangente unitaria
        t_xf = -distX / nLong
        t_yf = -distY / nLong
        # normal unitaria
        auxNormal_x.append(t_yf)
        auxNormal_y.append(-t_xf)
        
        iaux = iaux + 1

    aNormal_xf.append(auxNormal_x)                          # almacenar las normales
    aNormal_yf.append(auxNormal_y)

    aDist_cell_edge.append(aux_dist_cCell_cEdge)            # almacenar valores centro celda-lado

#obtenemos los maximos y los minimos de los centros de los bordes
nMax_Cx = np.amax(aCenint_x)
nMax_Cy = np.amax(aCenint_y)
nMin_Cx = np.amin(aCenint_x)
nMin_Cy = np.amin(aCenint_y)

# obtener los valores de distancia nodo - centro de celda
for vh in mesh.vertices():
    suma = 0
    aDistance = []
    vertex = mesh.point(vh)
    nVertex = vh.idx()
    for fh in mesh.vf(vh):
        nCell = fh.idx()
        disvX = aBarX[nCell] - vertex[0]
        disvY = aBarY[nCell] - vertex[1]
        dist = np.sqrt((disvX**2) + (disvY**2))
        aDistance.append(dist)
        aDist_node_cell[nVertex, nCell] = dist

    # calculo la suma de inversos de las distancias
    for valor in aDistance:
        suma = suma + pow(valor,-1)
    aSum_dist_node_cell.append(suma)                        # almaceno en el arreglo


############################################################################################################
############################################################################################################
# INFO CONECTIVIDAD
############################################################################################################

auxEdge = np.zeros(totEd)                                   # auxiliar que almacena si un lado es frontera (0,1)
aW_f = []                                                   # almacena datos para la interpolacion

for fh in mesh.faces():
    nCell = fh.idx()
    auxCell = []                                            # auxiliar que almacena las celdas vecinas
    fronteraux = []                                         # arreglo de bandera
    index = 0

    # Relacion entre las celda y sus lados
    for eh in mesh.fe(fh):
        nEdge = eh.idx()
        valueX = aCenint_x[nEdge]
        valueY = aCenint_y[nEdge]
        link_cell_to_edge[nCell, index] = nEdge

        # Evaluar la condidion de frontera
        if valueX == nMax_Cx or valueX == nMin_Cx or valueY == nMax_Cy or valueY == nMin_Cy:
            fronteraux.append(1)
            auxEdge[nEdge] = 1
        else:
            fronteraux.append(0)
            auxEdge[nEdge] = 0
        index = index + 1
    # print(nCell, fronteraux)

    # relleno el auxiliar con las caras vecinas
    for fh in mesh.ff(fh):
        nCell_v = fh.idx()
        auxCell.append(nCell_v)

    index = 0
    auxidx = 0                                              # recorre el auxiliar

    # Se establece la relacion incluyendo las fronteras
    for value in fronteraux:
        if value == 1:
            aux = int(-1)
            link_cell_to_cell[nCell, index] = aux
        else:
            link_cell_to_cell[nCell, index] = auxCell[auxidx]
            auxidx = auxidx + 1
        index = index + 1

    print(nCell, fronteraux)

#fin del for fh

# Calculo de W_f para la interpolacion en interfaces
for eh in mesh.edges():
    nEdge = eh.idx()
    dist_w = []                                             # guarda las distancias
    aux_edge_to_cell = []                                   # almacena las celdas adyacentes a la interface
    index = 0

    if auxEdge[nEdge] == 1:
        aux_edge_to_cell.append(-1)

    for cell in link_cell_to_edge:
        i = 0
        for edge in cell:
            # Calculo de la distacia para W_f
            if nEdge == edge:
                aux_edge_to_cell.append(index)
                dist_w.append(aDist_cell_edge[index][i])
            i = i + 1
        index = index + 1

    link_edge_to_cell.append(aux_edge_to_cell)              # almacenamos la relacion lado - celda

    # Calculo de W_f (Solo interiores)
    if auxEdge[nEdge] == 0:
        w_f = (1/dist_w[0]) / (1/dist_w[0] + 1/dist_w[1])
    else:
        w_f = 0                                             # No se asigna valor (frontera)

    aW_f.append(w_f)                                        # almaceno el valor

######################################## SOLO PARA TEST ####################################################

for fh in mesh.faces():
    nCell = fh.idx()
    aVt1 = []

    #guardo los vertices
    for vh in mesh.fv(fh):
        punto = mesh.point(vh)
        aVt1.append(punto)
    
    
    # gaurdo los puntos de interes
    for i in range(3):
        vecina = link_cell_to_cell[nCell][i]

        if vecina == -1:
            a = aVt1[i - 1]
            b = aVt1[i]

            #OK
            for fh in mesh.faces():
                nCell_p = fh.idx()
                aVt2 = []
                
                for vh in mesh.fv(fh):
                    punto1 = mesh.point(vh)
                    aVt2.append(punto1)

                for j in range(3):
                    vecina2 = link_cell_to_cell[nCell_p][j]
                    nEdge_p = link_cell_to_edge[nCell_p][j]

                    nEdge_p = int(nEdge_p)
                    vaY = aCenint_y[nEdge_p]
                    vaX = aCenint_x[nEdge_p]

                    if vecina2 == -1:
                        c = aVt2[j -1]
                        d = aVt2[j]

                        if c[0] == b[0] and d[0] == a[0]:
                            #print('Res: ', nCell_p,'c: ', c, 'd: ',d)
                            
                            if vaY == nMax_Cy or vaY == nMin_Cy:                                
                                celdaPar[nCell] = nCell_p
                                ladoPar[nCell] = nEdge_p

                        if c[1] == b[1] and d[1] == a[1]:

                            if vaX == nMax_Cx or vaX == nMin_Cx:                                
                                celdaPar[nCell] = nCell_p
                                ladoPar[nCell] = nEdge_p
    
    
        
print(celdaPar)       
  
############################################################################################################
############################################################################################################
# PARAMETROS
############################################################################################################

aU_vtc_PM = np.zeros(totVx)                                 # Concentracion de los vertices

sigma = 2
xc = -2
yc = -4
x0 = 2.5
y0 = 2.5
x1 = 0
y1 = 0
tau0 = 7.0
tau1 = 0.001
a0 = 1
cx = 1.5
cy = 1.5
lamb = 0.1
s = 1.0
epsilon = 0
alpha = 0.4
T0 = 0.2

coeff_D = 0.001
kappa = 0.01

tEnd = 3.0

# Concentracion inicial (celdas)
uPM = []
for fh in mesh.faces():
    initPM = fnc.exacta(aBarX[nCell], aBarY[nCell], 0, T0)
    uPM.append(initPM)

aUPM_timecell = uPM                                         # historico PM celdas
aUPM_timevertex = []                                        # historico PM vertices
aUPM_timeedge = []                                          # historico PM lados

hExacta = uPM                                               # historico exacta

aJ_B_10 = np.zeros(totEd)
aJ_B_25 = np.zeros(totEd)

########################################################################################################
# MANIPULACION DE DATOS

CFL = 0.5
#velc maxima
aVelcX = []
aVelcY = []

razon1 = []
aDT = []                            # de aqui escojo el dt
for fh in mesh.faces():
    nCell = fh.idx()
    vX = fnc.vx(aBarX[nCell],aBarY[nCell],0,T0)
    vY = fnc.vy(aBarX[nCell],aBarY[nCell],0,T0)

    valR = (np.sqrt(vX**2 + vY**2) + coeff_D) / aArea[nCell]
    razon1.append(valR)
    aVelcX.append(vX)
    aVelcY.append(vY)

for fh in mesh.faces():
    nCell = fh.idx()
    val = aArea[nCell] / max(razon1)
    aDT.append(val)

dt = CFL * min(aDT)                 # dt inicial


#arreglo de Tiempo
aTime = np.linspace(0,tEnd, 1200)
aTime = aTime.astype(np.int64)
# dt = 1 / len(aTime)


############################################################################################################
############################################################################################################
# SOLVER
############################################################################################################

hora = 0
dia = 1                                         # Lunes

t = 0

for val in aTime:
    u_Edges = np.zeros(totEd)                   # inicializo arreglo de concentraciones en interfaces
    aExacta = []                                # arreglo de valores para la exacta

    razon1 = []
    aDT = []

    ########################################################################################################
    # Calculo de los flujos vertice - nodo
    for vh in mesh.vertices():
        nVertex = vh.idx()
        aVxy = mesh.point(vh)
        #calculo el viento
        Vx = fnc.vx(aVxy[0],aVxy[1],t,T0)
        Vy = fnc.vy(aVxy[0],aVxy[1],t,T0)

        valR = (np.sqrt(Vx**2 + Vy**2) + coeff_D) / aArea[nCell]
        razon1.append(valR)

        aux_pm = []                             # auxiliar para sumar las concntraciones
        for fh in mesh.vf(vh):
            nCell = fh.idx()
            w_vi = (uPM[nCell] / aDist_node_cell[nVertex, nCell]) / aSum_dist_node_cell[nVertex]
            aux_pm.append(w_vi)

        aU_vtc_PM[nVertex] = sum(aux_pm)
    
    aUPM_timevertex.append(aU_vtc_PM)

    #######################################################################################################
    # Calculo de las concentraciones en las celdas
    for fh in mesh.faces():
        nCell = fh.idx()
        nQ_Edges = 3                            # todos son triangulos

        # almaceno la info de los vertices
        aVtx = []                               # x
        aVty = []                               # y
        aVtc = []                               # indice
        for vh in mesh.fv(fh):
            point = mesh.point(vh)
            aVtx.append(point[0])
            aVty.append(point[1])
            aVtc.append(vh.idx())

        # almacenamientos temporales
        aPM_diff = []
        aPM_adv = []

        # calculo flujos de interfaces
        for i in range(nQ_Edges):
            nEdge = int(link_cell_to_edge[nCell, i])                # interface
            nCell_v = int(link_cell_to_cell[nCell, i])              # celda vecina

            a = aVtc[i]
            b = aVtc[i-1]

            # concentracion en la interface
            u_Edges[nEdge] = aW_f[nEdge] * uPM[nCell] + ((1 - aW_f[nEdge]) * uPM[nCell_v])

            # flujo tangencial
            flux_tang_PM = (aU_vtc_PM[a] - aU_vtc_PM[b]) / aLongLados[nEdge]

            # G_f para flujo advectivo (Upwind)
            G_f = (Vx * aNormal_xf[nCell][i] + Vy * aNormal_yf[nCell][i]) * aLongLados[nEdge]

            if nCell_v == -1:
                # distancia delta_f
                nCell_p = int(celdaPar[nCell])
                nEdge_p = int(ladoPar[nCell])
                distx = np.abs(aBarX[nCell] - aCenint_x[nEdge])
                disty = np.abs(aBarY[nCell] - aCenint_y[nEdge])
                delta_distx = aNormal_xf[nCell][i] * distx
                delta_disty = aNormal_yf[nCell][i] * disty
                delta_dist = np.sqrt(pow(delta_distx,2) + pow(delta_disty,2))

                # calcular x_l que es el punto medio entre los centroides de las celdas
                #En el caso del Test tomo los puntos medio de la interface de frontera
                x_l = aCenint_x[nEdge]
                y_l = aCenint_y[nEdge]

                # t_f punto l_f (Flujo tangencial) Info GEOMETRICA
                tf_lf = (aVtx[i] - aVtx[i-1]) * (x_l - aBarX[nCell]) * (aVty[i] - aVty[i-1]) * (y_l - aBarY[nCell]) / aLongLados[nEdge]

                # valores flujos
                diff_f = coeff_D * (((fnc.exacta(aBarX[nCell],aBarY[nCell],t,T0) - uPM[nCell]) / delta_dist) - (flux_tang_PM * tf_lf / delta_dist)) * aLongLados[nEdge]
                # diff_f = coeff_D * ((fnc.exacta(aBarX[nCell],aBarY[nCell],t,T0) - uPM[nCell]) / delta_dist) * aLongLados[nEdge]
                adv_f = 0                                                                   # no influye
                
                aPM_diff.append(diff_f)
                #print('frontera', aPM_diff)
                aPM_adv.append(adv_f)
            else:
                # distancia delta_f
                distx = aBarX[nCell] - aBarX[nCell_v]
                disty = aBarY[nCell] - aBarY[nCell_v]
                delta_distx = aNormal_xf[nCell][i] * distx
                delta_disty = aNormal_yf[nCell][i] * disty
                delta_dist = np.sqrt(pow(delta_distx,2) + pow(delta_disty,2))

                # calcular x_l que es el punto medio entre los centroides de las celdas
                x_l = (aBarX[nCell] + aBarX[nCell_v]) / 2
                y_l = (aBarY[nCell] + aBarY[nCell_v]) / 2

                # t_f punto l_f (Flujo tangencial) Info GEOMETRICA
                tf_lf = (aVtx[i] - aVtx[i-1]) * (x_l - aBarX[nCell]) * (aVty[i] - aVty[i-1]) * (y_l - aBarY[nCell]) / aLongLados[nEdge]

                # Flujo Difusivo
                diff_f = coeff_D * (((uPM[nCell_v] - uPM[nCell]) / delta_dist) - (flux_tang_PM * tf_lf / delta_dist)) * aLongLados[nEdge]
                # diff_f = coeff_D * ((uPM[nCell_v] - uPM[nCell]) / delta_dist) * aLongLados[nEdge]

                #Flujo Advectivo
                adv_f = ((G_f + np.abs(G_f))/2 * uPM[nCell]) - ((G_f - np.abs(G_f))/2 * uPM[nCell_v])

                aPM_diff.append(diff_f)
                
                aPM_adv.append(adv_f)
            
            #### fin recorrer lados
        
        PM1 = uPM[nCell]                                          # valor anterior
        val_diff = sum(aPM_diff)                                  # difusion
        val_adv = sum(aPM_adv)                                    # adveccion
        srcPM = fnc.source(aBarX[nCell], aBarY[nCell], t, T0)     # emision

        #print('PM1',PM1)

        # Delta de Dirac
        d0x = fnc.delta(aBarX[nCell] - x0)
        d0y = fnc.delta(aBarY[nCell] - y0)
        d0t = fnc.delta(t - tau0)
        d1x = fnc.delta(aBarX[nCell] - x1)
        d1y = fnc.delta(aBarY[nCell] - y1)
        d1t = fnc.delta(t - tau1)

        area = aArea[nCell]

        # Concentracion de la celda
        PM2 = PM1 + dt * ((val_diff - val_adv) / area + srcPM) + dt * s *((d0x * d0y * d0t) + (d1x * d1y * d1t))
        uPM[nCell] = PM2
        valExacta = fnc.exacta(aBarX[nCell], aBarY[nCell], t, T0)
        aExacta.append(valExacta)
        
    # recalculo el dt
    for fh in mesh.faces():
        nCell = fh.idx()
        val = aArea[nCell] / max(razon1)
        aDT.append(val)

    dt = CFL * min(aDT)                 # dt para la prox iteracion
    
    #print('DT', dt, 'vX y Vy: ', Vx, Vy)
    
    ### paso por celda
    aUPM_timecell = np.vstack((aUPM_timecell, uPM))
    aUPM_timeedge.append(u_Edges)
    hExacta = np.vstack((hExacta, aExacta))

    t = t + dt

    hora = hora + 1

    #CICLO DE UN DIA
    if hora == 24:
        hora = 0
        dia = dia + 1

    # CICLO DE UNA SEMANA
    if dia == 8:
        dia = 1

    ### paso del tiempo


#####################################################################################################################################
# Exportar Datos
dfPM10 = pd.DataFrame(aUPM_timecell)
export = dfPM10.to_csv(r'PM_estimado.csv', index = None, header=True)

dfExac = pd.DataFrame(hExacta)
export = dfExac.to_csv(r'PM_exacto.csv', index = None, header=True)

    

print('Listo!')