
import numpy as np
import openmesh as om
import markov
import markovdir
import diffusive as difs
import viento
import time
import random
import pandas as pd
import interpolacion as intpoly
import ponderacion as pond

#### FUNCION PARA CAMBIAR COMA POR PUNTO DECIMAL
f = lambda x : (x.replace(",","."))

mesh = om.TriMesh()

# vertices de malla
vh0 = mesh.add_vertex([2061.98, 13970, 0])
vh1 = mesh.add_vertex([2061.98, 2602.69, 0])
vh2 = mesh.add_vertex([16028.02, 2602.69, 0])
vh3 = mesh.add_vertex([16028.02, 13970, 0])

vh4 = mesh.add_vertex([5120.30, 13970, 0])
vh5 = mesh.add_vertex([8871.72, 13970, 0])
vh6 = mesh.add_vertex([12995.69, 13970, 0])
vh7 = mesh.add_vertex([5120.30, 2602.69, 0])
vh8 = mesh.add_vertex([8871.72, 2602.69, 0])
vh9 = mesh.add_vertex([12995.69, 2602.69, 0])

vh10 = mesh.add_vertex([2061.98, 11356.19, 0])
vh11 = mesh.add_vertex([2061.98, 8241.86, 0])
vh12 = mesh.add_vertex([2061.98, 5694.78, 0])
vh13 = mesh.add_vertex([16028.02, 11356.19, 0])
vh14 = mesh.add_vertex([16028.02, 8241.86, 0])
vh15 = mesh.add_vertex([16028.02, 5694.78, 0])

vh16 = mesh.add_vertex([4661.12, 11945.68, 0])
vh17 = mesh.add_vertex([7858.06, 11789.97, 0])
vh18 = mesh.add_vertex([10873.06, 12068.03, 0])
vh19 = mesh.add_vertex([13368.23, 11667.62, 0])

vh20 = mesh.add_vertex([5631.47, 9487.59, 0])
vh21 = mesh.add_vertex([9183.62, 8931.46, 0])
vh22 = mesh.add_vertex([11531.51, 9676.67, 0])
vh23 = mesh.add_vertex([13463.53, 8564.41, 0])

vh24 = mesh.add_vertex([4123.97, 6039.58, 0])
vh25 = mesh.add_vertex([6983.02, 7196.33, 0])
vh26 = mesh.add_vertex([7615.47, 5272.12, 0])
vh27 = mesh.add_vertex([10751.77, 6139,68, 0])
vh28 = mesh.add_vertex([13385.56, 6017.33, 0])

vh29 = mesh.add_vertex([0, 13970, 0])
vh30 = mesh.add_vertex([0, 0, 0])
vh31 = mesh.add_vertex([18090, 0, 0])
vh32 = mesh.add_vertex([18090, 13970, 0])

vh33 = mesh.add_vertex([0, 10243.93, 0])
vh34 = mesh.add_vertex([0, 7062.86, 0])
vh35 = mesh.add_vertex([0, 4293.33, 0])
vh36 = mesh.add_vertex([18090, 10243.93, 0])
vh37 = mesh.add_vertex([18090, 7062.86, 0])
vh38 = mesh.add_vertex([18090, 4293.33, 0])

vh39 = mesh.add_vertex([2061.98, 0, 0])
vh40 = mesh.add_vertex([5120.30, 0, 0])
vh41 = mesh.add_vertex([8871.72, 0, 0])
vh42 = mesh.add_vertex([12995.69, 0, 0])
vh43 = mesh.add_vertex([16028.02, 0, 0])

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


# get all points of the mesh
point_array = mesh.points()

# write and read meshes
om.write_mesh('test.off', mesh)
mesh_2 = om.read_trimesh('test.off')

###### CALCULO DE CENTROIDES Y AREAS DE LAS CELDAS ######

#Recorrer todas las interfaces
contador = 0
for eh in mesh.edges():
    eIdx = eh.idx() # devuelve el indice global del elemento
    contador = contador + 1
print(contador)

# arreglo de baricentros
aCenterx = []
aCentery = []
# arreglo de Area
aArea = []
# arreglo de Longitud de interfaces (edges)
aLong = np.zeros(contador)
# arreglo de centros de las interfaces
aCenint_x = np.zeros(contador)
aCenint_y = np.zeros(contador)

# arreglo de normales a las interfaces
aNormal_xf = []
aNormal_yf = []

# dictancia nodo - centro de celda
aDist_node_cell = np.zeros((44,66))
aDist_cell_edge = []

# sumatoria de los inversos de las distancias nodo - centro (vertices)
aSum_dist_node_cell = []

# relacion celda a sus bordes
link_cell_to_edge = np.zeros((66,3))
link_edge_to_cell = []
# relacion celda a celda vecinas
link_cell_to_cell = np.zeros((66,3))
# relacion celda vertices
link_cell_to_vertex = np.zeros((66,3))

#############################################################################################################
#############################################################################################################

#iteramos sobre los vertices de la cara
for fh in mesh.faces():
    nCell = fh.idx()
    #inicializo el arreglo para los puntos de los nodos
    aVtx = []
    aVty = []
    index = 0
    for vh in mesh.fv(fh):
        point = mesh.point(vh)
        aVtx.append(point[0])
        aVty.append(point[1])
    # calculamos el baricentro de la celda
    sumax = aVtx[0] + aVtx[1] + aVtx[2]
    nCenterx = sumax / 3
    aCenterx.append(nCenterx)
    sumay = aVty[0] + aVty[1] + aVty[2]
    nCentery = sumay / 3
    aCentery.append(nCentery)
    # inicializo los vectores auxiliares para guardar la normal de la celda
    auxNormal_x = []
    auxNormal_y = []
    # calculamos el area de la celda
    nArea = aVtx[0] * (aVty[2] - aVty[1]) + aVtx[1] * (aVty[0] - aVty[2]) + aVtx[2] * (aVty[1] - aVty[0])
    nArea = 0.5 * np.absolute(nArea)
    aArea.append(nArea)

    #auxiliar para el calculo distancia centro celda - centro lado
    aux_dist_cCell_cEdge = []

    # recorremos los bordes de la celda
    for eh in mesh.fe(fh):
        nEdge = eh.idx()
        #calculamos el punto medio del borde
        nCenint_x = 0.5 * (aVtx[index - 1] + aVtx[index])
        nCenint_y = 0.5 * (aVty[index - 1] + aVty[index])
        aCenint_x[nEdge] = nCenint_x
        aCenint_y[nEdge] = nCenint_y
        #calculamos la longitud del segmento
        nLx = aVtx[index - 1] - aVtx[index]
        nLy = aVty[index - 1] - aVty[index]
        nLong = np.sqrt(pow(nLx,2) + pow(nLy,2))

        # Calcular la distancia baricentro a centro de segmento
        dif_x = aCenterx[nCell] - nCenint_x
        dif_y = aCentery[nCell] - nCenint_y
        distance = np.sqrt(pow(dif_x,2) + pow(dif_y,2))
        aux_dist_cCell_cEdge.append(distance)
        
        aLong[nEdge] = nLong
        index = index + 1
        # Obtenemos los vectores normales unitarios de cada interface (edge)
        # tangente unitaria
        t_xf = -nLx / nLong
        t_yf = -nLy / nLong
        # normal unitaria
        auxNormal_x.append(t_yf)
        auxNormal_y.append(-t_xf)
    aNormal_xf.append(auxNormal_x)
    aNormal_yf.append(auxNormal_y)

    #Las distancias centro a centro se guardan
    aDist_cell_edge.append(aux_dist_cCell_cEdge)

#obtenemos los maximos y los minimos de los centros de los bordes
nMax_Cx = np.amax(aCenint_x)
nMax_Cy = np.amax(aCenint_y)
nMin_Cx = np.amin(aCenint_x)
nMin_Cy = np.amin(aCenint_y)

print(aNormal_xf)
print(aNormal_yf)

# obtener los valores de distancia nodo - centro de celda
for vh in mesh.vertices():
    suma = 0
    aDistance = [] #auxiliar
    vertex = mesh.point(vh)
    vertex_x = vertex[0]
    vertex_y = vertex[1]
    nVertex = vh.idx()
    for fh in mesh.vf(vh):
        nCell = fh.idx()
        dif_x = aCenterx[nCell] - vertex_x
        dif_y = aCentery[nCell] - vertex_y
        # calculo la distancia
        distance = np.sqrt(pow(dif_x,2) + pow(dif_y,2))
        aDistance.append(distance)
        aDist_node_cell[nVertex, nCell] = distance

    for distancia in aDistance:
        suma = suma + pow(distancia,-1)
    #lo agrego a la sumatoria de inversos
    aSum_dist_node_cell.append(suma)

####INFO CONECTIVIDAD #######################################################
#############################################################################
idxEdge = np.zeros(contador)
for eh in mesh.edges():
    nEdge = eh.idx()
    idxEdge[nEdge] = nEdge
# auxiliar donde guardo si una interface corresponde a una frontera o no
auxEdge = np.zeros(contador)

# recorro las celdas 
for fh in mesh.faces():
    nCell = fh.idx()
    auxCell = []
    fronteraux = []
    index = 0
    for eh in mesh.fe(fh):
        nEdge = eh.idx()
        valueX = aCenint_x[nEdge]
        valueY = aCenint_y[nEdge]
        link_cell_to_edge[nCell, index] = nEdge

        if valueX == nMax_Cx or valueX == nMin_Cx or valueY == nMax_Cy or valueY == nMin_Cy:
            #evaluo la condicion de frontera
            fronteraux.append(1)
            auxEdge[nEdge] = 1
        else:
            fronteraux.append(0)
            auxEdge[nEdge] = 0
        index = index + 1

    for fh in mesh.ff(fh):
        nCell_v = fh.idx()
        auxCell.append(nCell_v)
    
    # OBTENEMOS EL VECTOR DE LINK DE LA CELDA CON SUS VECINOS
    index = 0
    auxidx = 0 # recorre el auxiliar
    for value in fronteraux:
        if value == 1:
            aux = int(-1)
            link_cell_to_cell[nCell, index] = aux
        else:
            link_cell_to_cell[nCell, index] = auxCell[auxidx]
            auxidx = auxidx + 1
        index = index + 1

    # OBTENEMOS LA MATRIZ DE LINK DE LA CELDA CON SUS VERTICES
    index = 0
    for vh in mesh.fv(fh):
        nVertex = vh.idx()
        link_cell_to_vertex[nCell,index] = nVertex
        index = index + 1

aW_f = []

for eh in mesh.edges():
    dist_w = [] # aqui agarro los valores de distancias correspondientes
    nEdge = eh.idx()
    # guarda las celdas adyacentes para una SOLA interface (auxiliar)
    aux_edge_to_cell = []
    index = 0
    if auxEdge[nEdge] == 1:
        aux_edge_to_cell.append(-1)
    
    for cell in link_cell_to_edge:
        i = 0
        for edge in cell:
            #  CALCULO DE W_f PARA INTERPOLACION DE INTERFACES
            if nEdge == edge:
                aux_edge_to_cell.append(index)
                dist_w.append(aDist_cell_edge[index][i])
            i = i + 1
        index = index + 1
    #finalmente guardamos los adyacentes en el global
    link_edge_to_cell.append(aux_edge_to_cell)

    # CALCULO W_F PARA :A INTERPOLACION. SOLO EN INTERIORES!!!
    if auxEdge[nEdge] == 0:
        w_f = (1/dist_w[0]) / (1/dist_w[0] + 1/dist_w[1])
    else:
        w_f = 0 # aun no asigno valor # se hace evaluando la condicion de borde
    
    aW_f.append(w_f)
    #print(dist_w)

#  CALCULO DE W_f PARA INTERPOLACION DE INTERFACES


####################################################################################################

####################################################################################################
## PARAMETROS
####################################################################################################

#concentraciones en vertices (lo inicializo como un arreglo de ceros)
aU_vtc_PM10 = np.zeros(44)
aU_vtc_PM25 = np.zeros(44)

# valores iniciales de concentracion (por ahora un numero random)
aUPM25_cell = np.array(random.sample(range(12,100), 66)) * pow(10,-9)
aUPM10_cell = np.array(random.sample(range(30,120), 66)) * pow(10,-9)
#aUPM25_cell = np.zeros(26)
#aUPM10_cell = np.zeros(26)


# arreglo de temperaturas
aTemp = []
# arreglo de velc viento
aVviento = []
# arreglo de dir viento
aDviento = []
# arreglo de lluvia
aLluvia = []

# emision residencial
srcRes_10 = np.array(random.sample(range(12,100), 66)) * pow(10,-9)
srcRes_25 = np.array(random.sample(range(12,100), 66)) * pow(10,-9)

# carga de emision por Fuentes Moviles en ruta (kg/h)
fMov_10 = np.array(random.sample(range(12,100), 66)) * pow(10,-9)
fMov_25 = np.array(random.sample(range(12,100), 66)) * pow(10,-9)

# carga de emision por Polvo Suspendido en ruta (kg/h)
pSusp_10 = np.array(random.sample(range(12,100), 66)) * pow(10,-9)
pSusp_25 = np.array(random.sample(range(12,100), 66)) * pow(10,-9)

# valores iniciales de concentracion (por ahora un numero random)
aUPM25_cell = []
aUPM10_cell = []

i = 0
for fh in mesh.faces():
    PM10val = srcRes_10[i]/(aArea[i]*200)
    PM25val = srcRes_25[i]/(aArea[i]*200)
    aUPM10_cell.append(PM10val)
    aUPM25_cell.append(PM25val)
    i = i + 1


# Matriz datos Valor/Tiempo
# PM 10
aUPM10_timecell = aUPM10_cell
aUPM10_timevertex = []
aUPM10_timeedge = []
# PM 2.5
aUPM25_timecell = aUPM25_cell
aUPM25_timevertex = []
aUPM25_timeedge = []

# Matrices J_B
aJ_B_10 = np.zeros(contador)
aJ_B_25 = np.zeros(contador)

####################################################################################################
## CARGA DE DATOS EN PANDA
####################################################################################################


####################################################################################################
# MANIPULACION DE DATOS

# the files with meteoriogical data are loaded in interpolacion file
temperatura = intpoly.resultado('temp')
lluvia = intpoly.resultado('rain')
vientoVel = intpoly.resultado('windVel')
vientoDir = intpoly.resultado('windDir')

largos = [np.size(temperatura),np.size(lluvia),np.size(vientoVel),np.size(vientoDir)]        
# Ahora defino el arreglo de tiempo
# tamanp de vector de tiempo
sizeT = np.size(temperatura) - 1
# Arreglo de tiempo
#aTime = np.linspace(0,sizeT,np.size(temperatura))
aTime = np.linspace(0,504, 505)
aTime = aTime.astype(np.int64)
dt = 1 / len(aTime)

totFlujA = []
totFlujD = []

####################################################################################################

####################################################################################################
####################################################################################################
#ESTO SE DEBE REPETIR EN Todo INTERVALO DE TIEMPO
####################################################################################################
hora = 0
dia = 1 #Lunes

for t in aTime:

    # Calculamos las componentes del viento
    V_x = viento.calc_Vx(vientoVel[t], vientoDir[t])
    V_y = viento.calc_Vy(vientoVel[t], vientoDir[t])

    # CALCULAR COEF DE DIFUSION

    #obtengo la lluvia
    nLluvia = lluvia[t]
    # obtengo la temp
    nTemp = temperatura[t] 
    # viscosidad del aire
    Nu = difs.visc_din(nTemp)
    # distancia de colision del aire
    Lmb = difs.freepath(nTemp)
    # Coef de difusion de PM 10 y 2.5 se lleva a horas
    Diff_cte1 = 2 * difs.coeffdif(nTemp,10.0,Lmb,Nu) * 3600
    Diff_cte2 = 2 * difs.coeffdif(nTemp,2.5,Lmb,Nu) * 3600

    # inicializo arreglo de concentraciones en las interfaces
    u25_Edges = np.zeros(contador)
    u10_Edges = np.zeros(contador)

    ##########################################
    # calcular el flujo en los vertices,nodo
    for vh in mesh.vertices():
        nVertex = vh.idx()
        #guardo la concentracion
        aux_pm25 = []
        aux_pm10 = []
        for fh in mesh.vf(vh):
            nCell = fh.idx()
            w_vi10 = (aUPM10_cell[nCell] / aDist_node_cell[nVertex, nCell]) / aSum_dist_node_cell[nVertex]
            w_vi25 = (aUPM25_cell[nCell] / aDist_node_cell[nVertex, nCell]) / aSum_dist_node_cell[nVertex]
            aux_pm10.append(w_vi10)
            aux_pm25.append(w_vi25)
        #la sumatoria y obtengo el valor de concentracion en cada vertice
        aU_vtc_PM10[nVertex] = sum(aux_pm10)
        aU_vtc_PM25[nVertex] = sum(aux_pm25)
        
    
    aUPM10_timevertex.append(aU_vtc_PM10)
    aUPM25_timevertex.append(aU_vtc_PM25)

    #para guardar los flujos en un arreglo (por ahora)
    flujosD = []
    flujosA = []

    #print('PM 25', aUPM25_cell)

    #######################################
    # Calcular la concentracion en cada celda
    for fh in mesh.faces():
        # inicializo la suma
        nCell = fh.idx()
        nQ_Edges = 3 # asumiendo q todos son triangulos
        ##inicializo el arreglo para los puntos de los nodos
        aVtx = []
        aVty = []
        aVtc = []
        index = 0

        # Matrices Aux PM 2.5
        aPM25_diff = []
        aPM25_adv = []
        # Matrices Aux PM 10
        aPM10_diff = []
        aPM10_adv = []

        for vh in mesh.fv(fh):
            point = mesh.point(vh)
            aVtx.append(point[0])
            aVty.append(point[1])
            aVtc.append(vh.idx())

        ####################################
        ## CALCULO FLUJOS DE LA INTERFACE
        for i in range(nQ_Edges):

            # obtenemos el borde de celda y la celda vecina
            nEdge = int(link_cell_to_edge[nCell, i])
            nCell_v = int(link_cell_to_cell[nCell, i])
            # obtengo los indices de los vertices
            a = aVtc[i]
            b = aVtc[i-1]

            # CALCULO LA CONCENTRACION EN LA INTERFACE
            u10_Edges[nEdge] = aW_f[nEdge] * aUPM10_cell[nCell] + (1 - aW_f[nEdge]) * aUPM10_cell[nCell_v]
            u25_Edges[nEdge] = aW_f[nEdge] * aUPM25_cell[nCell] + (1 - aW_f[nEdge]) * aUPM25_cell[nCell_v]

            #flujo tangencial (flujo vectores/L)
            flux_tang_PM10 = (aU_vtc_PM10[a] - aU_vtc_PM10[b]) / aLong[nEdge]
            flux_tang_PM25 = (aU_vtc_PM25[a] - aU_vtc_PM25[b]) / aLong[nEdge]

            # G_f para flujo advectivo (Upwind)
            G_f = (V_x * aNormal_xf[nCell][i] + V_y * aNormal_yf[nCell][i]) * aLong[nEdge]
            
            ###################################################################################
            ###################################################################################
            if nCell_v == -1:
                # distancia entre el centro de la celda y centro de la pared frontera
                distx = aCenterx[nCell] - aCenint_x[nEdge]
                disty = aCentery[nCell] - aCenint_y[nEdge]
                l_f = np.sqrt(pow(distx,2) + pow(disty,2))
                delta_distx = aNormal_xf[nCell][i] * distx
                delta_disty = aNormal_yf[nCell][i] * disty

                delta_dist = np.sqrt(pow(delta_distx,2) + pow(delta_disty,2))

                # calcular x_l que es el punto medio entre los centroides de las celdas
                x_l = (aCenterx[nCell] + aCenint_x[nEdge]) / 2
                y_l = (aCentery[nCell] + aCenint_y[nEdge]) / 2

                # t_f punto l_f (Flujo tangencial) Info GEOMETRICA
                tf_lf = (aVtx[i] - aVtx[i-1]) * (x_l - aCenterx[nCell]) * (aVty[i] - aVty[i-1]) * (y_l - aCentery[nCell]) / aLong[nEdge]

                # Para Primera iteracion
                if t == 0:
                    r_10 = aLong[nEdge] * Diff_cte1
                    
                    r_25 = aLong[nEdge] * Diff_cte2
                    
                    A = 0
                    aJ_B_10[nEdge] = r_10 * (A - aUPM10_cell[nCell])
                    aJ_B_25[nEdge] = r_25 * (A - aUPM25_cell[nCell])
                
                u10_Edges[nEdge] = aUPM10_cell[nCell] + aJ_B_10[nEdge] * delta_dist
                u25_Edges[nEdge] = aUPM25_cell[nCell] + aJ_B_25[nEdge] * delta_dist

                #recalculo J_B
                aJ_B_10[nEdge] = ((((u10_Edges[nEdge] - aUPM10_cell[nCell]) / delta_dist) - ((flux_tang_PM10 * tf_lf) / delta_dist))) * aLong[nEdge]
                aJ_B_25[nEdge] = ((((u25_Edges[nEdge] - aUPM25_cell[nCell]) / delta_dist) - ((flux_tang_PM25 * tf_lf) / delta_dist))) * aLong[nEdge]
                
                
                # FLUJO DIFUSIVO
                aux_diff_PM10 = aJ_B_10[nEdge]
                aux_diff_PM25 = aJ_B_25[nEdge]

                # FLUJO ADVECTIVO
                aux_adv_PM10 = 0
                aux_adv_PM25 = 0

                # AGREGO AMBOS FLUJOS A LOS ARREGLOS DE SUMATORIA
                aPM10_adv.append(aux_adv_PM10)
                aPM25_adv.append(aux_adv_PM25)
                aPM10_diff.append(aux_diff_PM10)
                aPM25_diff.append(aux_diff_PM25)
            else:
                # calculamos la distancia entre baricentros
                distx = aCenterx[nCell] - aCenterx[nCell_v]
                disty = aCentery[nCell] - aCentery[nCell_v]
                l_f = np.sqrt(pow(distx,2) + pow(disty,2))
                
                delta_distx = aNormal_xf[nCell][i] * distx
                delta_disty = aNormal_yf[nCell][i] * disty
                delta_dist = np.sqrt(pow(delta_distx,2) + pow(delta_disty,2))

                # calcular x_l que es el punto medio entre los centroides de las celdas
                x_l = (aCenterx[nCell] + aCenterx[nCell_v]) / 2
                y_l = (aCentery[nCell] + aCentery[nCell_v]) / 2

                # t_f punto l_f (segunda parte del flujo tangencial)
                tf_lf = (aVtx[i] - aVtx[i-1]) * (x_l - aCenterx[nCell]) * (aVty[i] - aVty[i-1]) * (y_l - aCentery[nCell]) / aLong[nEdge]

                # FLUJO DIFUSIVO
                aux_diff_PM10 = Diff_cte1 * ((((aUPM10_cell[nCell_v] - aUPM10_cell[nCell]) / delta_dist) - ((flux_tang_PM10 * tf_lf) / delta_dist))) * aLong[nEdge]
                aux_diff_PM25 = Diff_cte2 * ((((aUPM25_cell[nCell_v] - aUPM25_cell[nCell]) / delta_dist) - ((flux_tang_PM25 * tf_lf) / delta_dist))) * aLong[nEdge]

                # FLUJO ADVECTIVO
                aux_adv_PM10 = ((G_f + np.absolute(G_f)) * 0.5 * aUPM10_cell[nCell]) - ((G_f - np.absolute(G_f)) * 0.5 * aUPM10_cell[nCell_v])
                aux_adv_PM25 = ((G_f + np.absolute(G_f)) * 0.5 * aUPM25_cell[nCell]) - ((G_f - np.absolute(G_f)) * 0.5 * aUPM25_cell[nCell_v])

                # AGREGO AMBOS FLUJOS A LOS ARREGLOS DE SUMATORIA
                aPM10_adv.append(aux_adv_PM10)
                aPM25_adv.append(aux_adv_PM25)
                aPM10_diff.append(aux_diff_PM10)
                aPM25_diff.append(aux_diff_PM25)
        
            #FIN FLUJO INTERFACE

        # finalmente asigno el flujo a la celda
        PM10 = aUPM10_cell[nCell]
        PM25 = aUPM25_cell[nCell]
        value_diff10 = sum(aPM10_diff)
        value_diff25 = sum(aPM25_diff)
        value_adv10 = sum(aPM10_adv)
        value_adv25 = sum(aPM25_adv)

        flujosD.append(value_diff25)
        flujosA.append(value_adv25)


        # Velocidad de Deposicion
        vW25 = 4 * pow(10,4) * nLluvia * pow(10,-3)
        vW10 = 4 * pow(10,4) * nLluvia * pow(10,-3)

        area = aArea[nCell]

        # Emision
        src_PM10 = (srcRes_10[nCell] + fMov_10[nCell] + pSusp_10[nCell] * 0)/(area * 200)
        src_PM25 = (srcRes_25[nCell] + fMov_25[nCell] + pSusp_25[nCell] * 0)/(area * 200)

        aUPM10_cell[nCell] = PM10 + dt * (value_diff10 - value_adv10) / area - dt * (vW10 * PM10) + dt * src_PM10 * pond.calcular(hora, 10)
        aUPM25_cell[nCell] = PM25 + dt * (value_diff25 - value_adv25) / area - dt * (vW25 * PM25) + dt * src_PM25 * pond.calcular(hora, 2)

        if aUPM10_cell[nCell] < 0:
            aUPM10_cell[nCell] = 0

        if aUPM25_cell[nCell] < 0:
            aUPM25_cell[nCell] = 0

        #print('Difusion: ', value_diff10)
        #print('Adveccion', value_adv10)
        
        # FIN FLUJO EN TODAS LAS CARAS
        
    # print('HOLAAA', aUPM10_cell)

    totFlujA.append(flujosA)
    totFlujD.append(flujosD)
    b = aUPM10_cell
    aUPM10_timecell = np.vstack((aUPM10_timecell,aUPM10_cell))
    aUPM25_timecell = np.vstack((aUPM25_timecell,aUPM25_cell))

    aUPM10_timeedge.append(u10_Edges)
    aUPM25_timeedge.append(u25_Edges)

    hora = hora + 1

    #CICLO DE UN DIA
    if hora == 24:
        hora = 0
        dia = dia + 1

    # CICLO DE UNA SEMANA
    if dia == 8:
        dia = 1
    


#print('Viento',aVviento)

##### CONVERIR TOdO A M2!!!

dfPM25 = pd.DataFrame(aUPM25_timecell)
dfPM10 = pd.DataFrame(aUPM10_timecell)
temperatura = pd.DataFrame(aTemp)
dfViento = pd.DataFrame(aVviento)
dfLluvia = pd.DataFrame(aLluvia)

export = dfPM25.to_csv(r'PM25_data.csv', index = None, header=True)
export = dfPM10.to_csv(r'PM10_data.csv', index = None, header=True)

#######################################################################
# dfAdv10 = pd.DataFrame(totFlujA)
# dfDif10 = pd.DataFrame(totFlujD)

# export = dfAdv10.to_csv(r'PM25_Adv.csv', index = None, header=True)
# export = dfDif10.to_csv(r'PM25_Dif.csv', index = None, header=True)

print('OK')

#ESTABILIDAD
est = max(aArea)/(max(vientoVel) * 3600)
print('Estabilidad: ',est)
print('Delta t: ',len(aTime))