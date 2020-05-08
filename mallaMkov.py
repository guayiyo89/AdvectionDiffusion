import numpy as np
import openmesh as om
import markov
import markovdir
import diffusive as difs
import viento
import time
import random
import pandas as pd

#### FUNCION PARA CAMBIAR COMA POR PUNTO DECIMAL
f = lambda x : (x.replace(",","."))


mesh = om.TriMesh()

# vertices de malla
vh0 = mesh.add_vertex([0, 11356.46, 0])
vh1 = mesh.add_vertex([6479, 11356.46, 0])
vh2 = mesh.add_vertex([12254, 11356.46, 0])
vh3 = mesh.add_vertex([17732, 11356.46, 0])
vh4 = mesh.add_vertex([17732, 4744.824, 0])
vh5 = mesh.add_vertex([17732, 0, 0])
vh6 = mesh.add_vertex([8712, 0, 0])
vh7 = mesh.add_vertex([0, 0, 0])
vh8 = mesh.add_vertex([0, 4900.392, 0])
vh9 = mesh.add_vertex([0, 7467.264, 0])
# vertices interiores
vh10 = mesh.add_vertex([3905, 9022.944, 0])
vh11 = mesh.add_vertex([8129, 8567.352, 0])
vh12 = mesh.add_vertex([5808, 5678.232, 0])
vh13 = mesh.add_vertex([10692, 6311.616, 0])
vh14 = mesh.add_vertex([13607, 8711.808, 0])
vh15 = mesh.add_vertex([15136, 6900.552, 0])
vh16 = mesh.add_vertex([3553, 2833.56, 0])
vh17 = mesh.add_vertex([7568, 2922.456, 0])
vh18 = mesh.add_vertex([13024, 3411.384, 0])

# add a couple of faces to the mesh
fh0 = mesh.add_face(vh0, vh10, vh1)
fh1 = mesh.add_face(vh0, vh9, vh10)
fh2 = mesh.add_face(vh10, vh9, vh12)
fh3 = mesh.add_face(vh12, vh9, vh8)
fh4 = mesh.add_face(vh8, vh16, vh12)
fh5 = mesh.add_face(vh8, vh7, vh16)
fh6 = mesh.add_face(vh7, vh6, vh16)
fh7 = mesh.add_face(vh12, vh16, vh17)
fh8 = mesh.add_face(vh12, vh11, vh10)
fh9 = mesh.add_face(vh11, vh12, vh13)
fh10 = mesh.add_face(vh13, vh12, vh17)
fh11 = mesh.add_face(vh13, vh17, vh18)
fh12 = mesh.add_face(vh13, vh18, vh15)
fh13 = mesh.add_face(vh13, vh15, vh14)
fh14 = mesh.add_face(vh3, vh14, vh15)
fh15 = mesh.add_face(vh14, vh11, vh13)
fh16 = mesh.add_face(vh1, vh10, vh11)
fh17 = mesh.add_face(vh1, vh11, vh2)
fh18 = mesh.add_face(vh2, vh11, vh14)
fh19 = mesh.add_face(vh2, vh14, vh3)
fh20 = mesh.add_face(vh3, vh15, vh4)
fh21 = mesh.add_face(vh4, vh15, vh18)
fh22 = mesh.add_face(vh4, vh18, vh5)
fh23 = mesh.add_face(vh17, vh16, vh6)
fh24 = mesh.add_face(vh18, vh17, vh6)


vh_list = [vh5, vh18, vh6]
fh25 = mesh.add_face(vh_list)

# get all points of the mesh
point_array = mesh.points()

# write and read meshes
om.write_mesh('test.off', mesh)
mesh_2 = om.read_trimesh('test.off')

###### CALCULO DE CENTROIDES Y AREAS DE LAS CELDAS ######

# arreglo de baricentros
aCenterx = []
aCentery = []
# arreglo de Area
aArea = []
# arreglo de Longitud de interfaces (edges)
aLong = np.zeros(44)
# arreglo de centros de las interfaces
aCenint_x = np.zeros(44)
aCenint_y = np.zeros(44)

# arreglo de normales a las interfaces
aNormal_xf = []
aNormal_yf = []

# dictancia nodo - centro de celda
aDist_node_cell = np.zeros((19,26))

# sumatoria de los inversos de las distancias nodo - centro (vertices)
aSum_dist_node_cell = []

# relacion celda a sus bordes
link_cell_to_edge = np.zeros((26,3))
# relacion celda a celda vecinas
link_cell_to_cell = np.zeros((26,3))
# relacion celda vertices
link_cell_to_vertex = np.zeros((26,3))

#############################################################################################################
#############################################################################################################

#iteramos sobre los vertices de la cara
for fh in mesh.faces():
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
    nCenterx = 0.5 * sumax
    aCenterx.append(nCenterx)
    sumay = aVty[0] + aVty[1] + aVty[2]
    nCentery = 0.5 * sumay
    aCentery.append(nCentery)
    # inicializo los vectores auxiliares para guardar la normal de la celda
    auxNormal_x = []
    auxNormal_y = []
    # calculamos el area de la celda
    nArea = aVtx[0] * (aVty[2] - aVty[1]) + aVtx[1] * (aVty[0] - aVty[2]) + aVtx[2] * (aVty[1] - aVty[0])
    nArea = 0.5 * np.absolute(nArea)
    aArea.append(nArea)
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
        aLong[nEdge] = nLong
        index = index + 1
        # Obtenemos los vectores normales unitarios de cada interface (edge)
        # tangente unitaria
        t_xf = -nLx / nLong
        t_yf = -nLy / nLong
        # normal unitaria (cambie los signos pues la malla recorre al reves)
        auxNormal_x.append(t_xf)
        auxNormal_y.append(-t_xf)
    aNormal_xf.append(auxNormal_x)
    aNormal_yf.append(auxNormal_y)

#obtenemos los maximos y los minimos de los centros de los bordes
nMax_Cx = np.amax(aCenint_x)
nMax_Cy = np.amax(aCenint_y)
nMin_Cx = np.amin(aCenint_x)
nMin_Cy = np.amin(aCenint_y)


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
        else:
            fronteraux.append(0)
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

####################################################################################################
####################################################################################################
## PARAMETROS
####################################################################################################

#concentraciones en vertices (lo inicializo como un arreglo de ceros)
aU_vtc_PM10 = np.zeros((19,26))
aU_vtc_PM25 = np.zeros((19,26))

# VIENTO (m/h) (corrdenadas x e y) (HACER FUNCION QUE DESCOMPONGA ESTE VALOR)
velc = 4500
dirviento = 230.0

# Precipitaciones inicial
nLluvia = 0.98

# Temperatura Inicial (Promedio)
nTemp = 11.0

# Arreglo de tiempo
aTime = np.linspace(0,768,768)

# valores iniciales de concentracion (por ahora un numero random)
aUPM25_cell = np.array(random.sample(range(12,100), 26)) * pow(10,-9)
aUPM10_cell = np.array(random.sample(range(30,120), 26)) * pow(10,-9)

# valores obtenido de las medianas de los datos de PM (para establecer condicion de emision)
aUPM25_med = np.array(random.sample(range(12,55), 26)) * pow(10,-9)
aUPM10_med = np.array(random.sample(range(25,70), 26)) * pow(10,-9)

# Matriz datos Valor/Tiempo
# PM 10
aUPM10_timecell = aUPM10_cell
aUPM10_timevertex = []
# PM 2.5
aUPM25_timecell = aUPM25_cell
print(aUPM25_timecell)
aUPM25_timevertex = []

# arreglo de temperaturas
aTemp = [nTemp]
# arreglo de velc viento
aVviento = [velc]
# arreglo de dir viento
aDviento = []
# arreglo de lluvia
aLluvia = []

dt = 1 / len(aTime)


####################################################################################################
####################################################################################################
#ESTO SE DEBE REPETIR EN Todo INTERVALO DE TIEMPO
####################################################################################################


for time in aTime:
    # Calculamos las componentes del viento
    V_x = viento.calc_Vx(velc, dirviento)
    V_y = viento.calc_Vy(velc, dirviento)

    # CALCULAR COEF DE DIFUSION
    # viscosidad del aire
    Nu = difs.visc_din(nTemp)
    # distancia de colision del aire
    Lmb = difs.freepath(nTemp)
    # Coef de difusion de PM 10 y 2.5 se lleva a horas
    Diff_cte1 = difs.coeffdif(nTemp,10.0,Lmb,Nu) * 3600
    Diff_cte2 = difs.coeffdif(nTemp,2.5,Lmb,Nu) * 3600

    ##########################################
    # calcular el flujo en los vertices,nodo
    for vh in mesh.vertices():
        nVertex = vh.idx()
        for fh in mesh.vf(vh):
            nCell = fh.idx()
            valor_u10 = (aUPM10_cell[nCell] / aDist_node_cell[nVertex, nCell]) / aSum_dist_node_cell[nVertex]
            aU_vtc_PM10[nVertex, nCell] = valor_u10
            valor_u25 = (aUPM25_cell[nCell] / aDist_node_cell[nVertex, nCell]) / aSum_dist_node_cell[nVertex]
            aU_vtc_PM25[nVertex, nCell] = valor_u25
    
    aUPM10_timevertex.append(aU_vtc_PM10)
    aUPM25_timevertex.append(aU_vtc_PM25)

    #print('PM 25', aUPM25_cell)

    #######################################
    # Calcular el flujo en cada celda
    for fh in mesh.faces():
        # inicializo la suma
        nUPM10_cell = 0 # suma para el flujo difusivo
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
            # obtengo los vertices
            a = aVtc[i - 1]
            b = aVtc[i]

            #flujo tangencial
            flux_tang_PM10 = (aU_vtc_PM10[a, nCell] - aU_vtc_PM10[b, nCell]) / aLong[nEdge]
            flux_tang_PM25 = (aU_vtc_PM25[a, nCell] - aU_vtc_PM25[b, nCell]) / aLong[nEdge]

            # G_f para flujo advectivo (Upwind)
            G_f = (V_x * aNormal_xf[nCell][i] + V_y * aNormal_yf[nCell][i]) * aLong[nEdge]
            
            ###################################################################################
            ###################################################################################
            if nCell_v == -1:
                # distancia entre el centro de la celda y centro de la pared frontera
                distx = aCenterx[nCell] - aCenint_x[nEdge]
                disty = aCentery[nCell] - aCenint_y[nEdge]
                delta_dist = np.sqrt(pow(distx,2) + pow(disty,2))

                # calcular x_l que es el punto medio entre los centroides de las celdas
                x_l = (aCenterx[nCell] + aCenint_x[nEdge]) / 2
                y_l = (aCentery[nCell] + aCenint_y[nEdge]) / 2

                # t_f punto l_f (segunda parte del flujo tangencial) Info GEOMETRICA
                tf_lf = (aVtx[i-1] - aVtx[i]) * (x_l - aCenterx[nCell]) * (aVty[i-1] - aVty[i]) * (y_l - aCentery[nCell]) / aLong[nEdge]
                
                # FLUJO DIFUSIVO
                aux_diff_PM10 = Diff_cte1 * (((0 - aUPM10_cell[nCell]) / delta_dist) - ((flux_tang_PM10 * tf_lf) / delta_dist)) * aLong[nEdge]
                aux_diff_PM25 = Diff_cte2 * (((0 - aUPM25_cell[nCell]) / delta_dist) - ((flux_tang_PM25 * tf_lf) / delta_dist)) * aLong[nEdge]

                # FLUJO ADVECTIVO
                aux_adv_PM10 = (G_f + np.absolute(G_f)) * 0.5 * aUPM10_cell[nCell]
                aux_adv_PM25 = (G_f + np.absolute(G_f)) * 0.5 * aUPM25_cell[nCell]

                # AGREGO AMBOS FLUJOS A LOS ARREGLOS DE SUMATORIA
                aPM10_adv.append(aux_adv_PM10)
                aPM25_adv.append(aux_adv_PM25)
                aPM10_diff.append(aux_diff_PM10)
                aPM25_diff.append(aux_diff_PM25)
            else:
                # calculamos la distancia entre baricentros
                distx = aCenterx[nCell] - aCenterx[nCell_v]
                disty = aCentery[nCell] - aCentery[nCell_v]
                delta_dist = np.sqrt(pow(distx,2) + pow(disty,2))

                # calcular x_l que es el punto medio entre los centroides de las celdas
                x_l = (aCenterx[nCell] + aCenterx[nCell_v]) / 2
                y_l = (aCentery[nCell] + aCentery[nCell_v]) / 2

                # t_f punto l_f (segunda parte del flujo tangencial)
                tf_lf = (aVtx[i-1] - aVtx[i]) * (x_l - aCenterx[nCell]) * (aVty[i-1] - aVty[i]) * (y_l - aCentery[nCell]) / aLong[nEdge]

                # FLUJO DIFUSIVO
                aux_diff_PM10 = Diff_cte1 * ((((aUPM10_cell[nCell_v] - aUPM10_cell[nCell]) / delta_dist) - ((flux_tang_PM10 * tf_lf) / delta_dist)) * aLong[nEdge]) * aLong[nEdge]
                aux_diff_PM25 = Diff_cte2 * ((((aUPM25_cell[nCell_v] - aUPM25_cell[nCell]) / delta_dist) - ((flux_tang_PM25 * tf_lf) / delta_dist)) * aLong[nEdge]) * aLong[nEdge]

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

        # Deposicion Seca (vels en m/s lelvarlo a m/h)
        dseca25 = difs.vdep(nTemp,2.5,aUPM25_cell[nCell],velc) * aUPM25_cell[nCell] * 900
        dseca10 = difs.vdep(nTemp,10,aUPM10_cell[nCell],velc) * aUPM10_cell[nCell] * 900

        # Emision
        cfSrc = difs.pond(nTemp) * 1.1

        aUPM10_cell[nCell] = PM10 + dt * (value_diff10 - value_adv10) / aArea[nCell] - dt * (nLluvia * pow(10,-3) * aUPM10_cell[nCell]/4 + dseca10) + cfSrc * aUPM10_med[nCell]
        aUPM25_cell[nCell] = PM25 + dt * (value_diff25 - value_adv25) / aArea[nCell] - dt * (nLluvia * pow(10,-3) * aUPM25_cell[nCell]/4 + dseca25) + cfSrc * aUPM25_med[nCell]
        # FIN FLUJO EN TODAS LAS CARAS
    # print('HOLAAA', aUPM10_cell)
    b = aUPM10_cell
    aUPM10_timecell = np.vstack((aUPM10_timecell,aUPM10_cell))
    aUPM25_timecell = np.vstack((aUPM25_timecell,aUPM25_cell))

    # Actualizo valores de Temperatura, Velocidad y Dir Viento
    nTemp = markov.temperatura(nTemp)

    velc = markov.valor_velc(velc/1000)
    # lo recnvierto a escala el viento
    velc = velc * 1000

    dirviento = markovdir.direction(dirviento)
    nLluvia = markovdir.ppm(nLluvia)

    # guardo valores
    aTemp.append(nTemp)
    aVviento.append(velc)
    aDviento.append(dirviento)
    aLluvia.append(nLluvia)


#print('Viento',aVviento)

##### CONVERIR TOdO A M2!!!

print(dseca25)

dfPM25 = pd.DataFrame(aUPM25_timecell)
dfPM10 = pd.DataFrame(aUPM10_timecell)
temperatura = pd.DataFrame(aTemp)
dfViento = pd.DataFrame(aVviento)
dfLluvia = pd.DataFrame(aLluvia)

print('PM 25', dfPM25)

export = dfPM25.to_csv(r'PM25_test.csv', index = None, header=True)
export = dfPM10.to_csv(r'PM10_test.csv', index = None, header=True)
export = temperatura.to_csv(r'temperatura_test.csv', index = None, header=True)
export = dfLluvia.to_csv(r'lluvia_test.csv', index = None, header=True)
export = dfViento.to_csv(r'viento_test.csv', index = None, header=True)