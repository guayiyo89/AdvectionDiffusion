import numpy as np
import openmesh as om
import pandas as pd
import exactas as fnc
from numba import jit

#### Defino el tipo de Malla
mesh = om.TriMesh() # TriMesh es una funcion de OpenMesh (om)

mesh = om.read_trimesh('prueba1000cell_2.off')
#mesh = om.read_trimesh('dmz.off')

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

link_cell_to_vertex = []                                # relacion celda vertices(indices)
link_cell_to_VtX = []
link_cell_to_VtY = []

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
        link_cell_to_edge[nCell, index] = int(nEdge)

        #evaluar condicion de frontera
        if mesh.is_boundary(eh):
            fronteraux.append(1)
            auxEdge[nEdge] = 1
        else:
            fronteraux.append(0)
            auxEdge[nEdge] = 0
        index += 1

    # relleno el auxiliar con las caras vecinas
    for fh in mesh.ff(fh):
        nCell_v = fh.idx()
        auxCell.append(int(nCell_v))

    index = 0
    auxidx = 0                                              # recorre el auxiliar

    # Se establece la relacion incluyendo las fronteras
    for value in fronteraux:
        if value == 1:
            link_cell_to_cell[nCell, index] = int(-1)
        else:
            link_cell_to_cell[nCell, index] = auxCell[auxidx]
            auxidx = auxidx + 1
        index = index + 1

    #print(nCell, fronteraux)

# obtencion de los vertices de cada una de las caras
for fh in mesh.faces():
    # almaceno la info de los vertices
    aVtx = []                               # x
    aVty = []                               # y
    aVtc = []                               # indice
    for vh in mesh.fv(fh):
        point = mesh.point(vh)
        aVtx.append(point[0])
        aVty.append(point[1])
        aVtc.append(vh.idx())
    link_cell_to_VtX.append(aVtx)
    link_cell_to_VtY.append(aVty)
    link_cell_to_vertex.append(aVtc)


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

nMax_Cx = np.amax(aCenint_x)
nMax_Cy = np.amax(aCenint_y)
nMin_Cx = np.amin(aCenint_x)
nMin_Cy = np.amin(aCenint_y)

for fh in mesh.faces():
    nCell = fh.idx()
    aVt1 = []

    if mesh.is_boundary(fh):

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

                    if mesh.is_boundary(fh):
                    
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
#tau0 = 7.0
#tau1 = 0.22
tau0 = 0.08
tau1 = 0.18
a0 = 1
cx = 1.5
cy = 1.5
lamb = 0.1
s = 0
epsilon = 0
alpha = 0.4
T0 = 0.2

coeff_D = 0.01
kappa = 0.01

tEnd = 3.0

# Concentracion inicial (celdas)
uPM = []
for fh in mesh.faces():
    nCell = fh.idx()
    #initPM = fnc.exacta(aBarX[nCell], aBarY[nCell], 0, T0, x0, x1, y0, y1, tau0, tau1)
    initPM = fnc.exacta2(aBarX[nCell], aBarY[nCell], 0, T0, x0, x1, y0, y1, tau0, tau1, 0)
    uPM.append(initPM)

aUPM_timecell = uPM                                         # historico PM celdas
aUPM_timevertex = []                                        # historico PM vertices
aUPM_timeedge = []                                          # historico PM lados

hExacta = uPM                                               # historico exacta

########################################################################################################
# MANIPULACION DE DATOS

CFL = 1
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


############################################################################################################
############################################################################################################
# SOLVER
############################################################################################################
#Funciones

@jit
def dist_delta(nCell,i, distx, disty):
    delta_distx = aNormal_xf[nCell][i] * distx
    delta_disty = aNormal_yf[nCell][i] * disty
    delta_dist = np.sqrt(pow(delta_distx,2) + pow(delta_disty,2))
    return delta_dist

@jit
def flux_tg(nCell, nEdge, i, a, b, x_l, y_l):
    #busqueda de la posicion de los vertices
    vetX_a = link_cell_to_VtX[nCell][i]
    vetY_a = link_cell_to_VtY[nCell][i]
    vetX_b = link_cell_to_VtX[nCell][i-1]
    vetY_b = link_cell_to_VtY[nCell][i-1]

    flux_tang_PM = (aU_vtc_PM[a] - aU_vtc_PM[b]) / aLongLados[nEdge]
    tf_lf = (vetX_a - vetX_b) * (x_l - aBarX[nCell]) * (vetY_a - vetY_b) * (y_l - aBarY[nCell]) / aLongLados[nEdge]
    flux_tan = flux_tang_PM * tf_lf
    return flux_tan

@jit
def gF(nCell, nEdge, i, Vx, Vy):
    g_f = (Vx * aNormal_xf[nCell][i] + Vy * aNormal_yf[nCell][i]) * aLongLados[nEdge]
    return g_f

@jit
def adveccion(G_f, nCell_v, nCell):
    adv = ((G_f + np.abs(G_f))/2 * uPM[nCell]) + ((G_f - np.abs(G_f))/2 * uPM[nCell_v])
    return adv

@jit
def difusion(nCell, nCell_v, nEdge, tf_lf, delta_dist):
    diff = coeff_D * (((uPM[nCell_v] - uPM[nCell]) / delta_dist) - (tf_lf / delta_dist)) * aLongLados[nEdge]
    return diff

@jit
def f_semidiscreto(nCell, uPM, aUPM_vtc, t):
    global aW_f

    #calculo el viento
    Vx = fnc.vx(aBarX[nCell],aBarY[nCell],t,T0)
    Vy = fnc.vy(aBarX[nCell],aBarY[nCell],t,T0)

    # almacenamientos temporales
    aPM_diff = []
    aPM_adv = []

    #calculo flujo de interfaces
    for i in range(3):
        nEdge = int(link_cell_to_edge[nCell, i])                # interface
        nCell_v = int(link_cell_to_cell[nCell, i])              # celda vecina

        a = link_cell_to_vertex[nCell][i]
        b = link_cell_to_vertex[nCell][i-1]

        # concentracion en la interface
        u_Edges[nEdge] = aW_f[nEdge] * uPM[nCell] + ((1 - aW_f[nEdge]) * uPM[nCell_v])

        if nCell_v == -1:
            nCell_p = int(celdaPar[nCell])                      # celda PAR
            nEdge_p = int(ladoPar[nCell])
            distx = np.abs(aBarX[nCell] - aCenint_x[nEdge]) + np.abs(aBarX[nCell_p] - aCenint_x[nEdge_p])
            disty = np.abs(aBarY[nCell] - aCenint_y[nEdge]) + np.abs(aBarY[nCell_p] - aCenint_y[nEdge_p])
            delta_dis = dist_delta(nCell, i, distx, disty)

            #En el caso del Test tomo los puntos medio de la interface de frontera
            x_l = aCenint_x[nEdge]
            y_l = aCenint_y[nEdge]

            t_f = flux_tg(nCell, nEdge, i, a, b, x_l, y_l)

            #flujo difusivo
            diff_f = difusion(nCell,nCell_p, nEdge, t_f, delta_dis)
            #flujo advectivo
            G_f = gF(nCell, nEdge, i, Vx, Vy)
            adv_f = adveccion(G_f,nCell_p,nCell)

            aPM_diff.append(diff_f)
            aPM_adv.append(adv_f)
        else:
            distx = aBarX[nCell] - aBarX[nCell_v]
            disty = aBarY[nCell] - aBarY[nCell_v]
            delta_dis = dist_delta(nCell, i, distx, disty)

            # calcular x_l que es el punto medio entre los centroides de las celdas
            x_l = (aBarX[nCell] + aBarX[nCell_v]) / 2
            y_l = (aBarY[nCell] + aBarY[nCell_v]) / 2

            t_f = flux_tg(nCell, nEdge, i, a, b, x_l, y_l)
            #flujo difusivo
            diff_f = difusion(nCell,nCell_v, nEdge, t_f, delta_dis)
            #flujo advectivo
            G_f = gF(nCell, nEdge, i, Vx, Vy)
            adv_f = adveccion(G_f,nCell_v,nCell)

            aPM_diff.append(diff_f)
            aPM_adv.append(adv_f)

        #### fin recorrer lados
    val_diff = sum(aPM_diff)                                  # difusion
    val_adv = sum(aPM_adv)                                    # adveccion

    #srcPM = fnc.source(aBarX[nCell], aBarY[nCell], t, T0, x0, x1, y0, y1, tau0, tau1, kappa)
    srcPM = fnc.source2(aBarX[nCell], aBarY[nCell], t, T0, x0, x1, y0, y1, tau0, tau1, kappa, hora)

    # Delta de Dirac
    d0x = fnc.delta(aBarX[nCell] - x0)
    d0y = fnc.delta(aBarY[nCell] - y0)
    #d0t = fnc.delta(t - tau0)
    d0t = fnc.delta(hora - tau0)
    d1x = fnc.delta(aBarX[nCell] - x1)
    d1y = fnc.delta(aBarY[nCell] - y1)
    #d1t = fnc.delta(t - tau1) 
    d1t = fnc.delta(hora - tau1) 

    res = (( (val_diff - val_adv) / aArea[nCell]) + srcPM) + s *((d0x * d0y * d0t) + (d1x * d1y * d1t)) - kappa * uPM[nCell] #sumae +srcPM
    return res
    
#############################################################################################################
            
t = 0
hora = 0
aTime = []                #guardo los val de tiempo
while t <= tEnd:

    u_Edges = np.zeros(totEd)                   # inicializo arreglo de concentraciones en interfaces
    aExacta = []                                # arreglo de valores para la exacta

    razon1 = []
    aDT = []

    ########################################################################################################
    # Calculo de los flujos vertice - nodo
    for vh in mesh.vertices():
        nVertex = vh.idx()
        aVxy = mesh.point(vh)

        aux_pm = []                             # auxiliar para sumar las concntraciones
        for fh in mesh.vf(vh):
            nCell = fh.idx()
            w_vi = (uPM[nCell] / aDist_node_cell[nVertex, nCell]) / aSum_dist_node_cell[nVertex]
            aux_pm.append(w_vi)

        aU_vtc_PM[nVertex] = sum(aux_pm)
    
    aUPM_timevertex.append(aU_vtc_PM)

    #######################################################################################################
    # Calculo de las concentraciones en las celdas

    #solver
    for fh in mesh.faces():
        nCell = fh.idx()
        PM1 = uPM[nCell]
        PMf = f_semidiscreto(nCell, uPM, aU_vtc_PM, t)
        PM2 = PM1 + dt * PMf
        uPM[nCell] = PM2
        #valExacta = fnc.exacta(aBarX[nCell], aBarY[nCell], t, T0, x0, x1, y0, y1, tau0, tau1)
        valExacta = fnc.exacta2(aBarX[nCell], aBarY[nCell], t, T0, x0, x1, y0, y1, tau0, tau1, hora)
        aExacta.append(valExacta)

    ### paso por celda
    aUPM_timecell = np.vstack((aUPM_timecell, uPM))
    aUPM_timeedge.append(u_Edges)
    hExacta = np.vstack((hExacta, aExacta))

    # recalculo el dt
    for fh in mesh.faces():
        nCell = fh.idx()
        vX = fnc.vx(aBarX[nCell],aBarY[nCell],t,T0)
        vY = fnc.vy(aBarX[nCell],aBarY[nCell],t,T0)
        valR = (np.sqrt(vX**2 + vY**2) + coeff_D) / aArea[nCell]
        razon1.append(valR)

    for fh in mesh.faces():
        nCell = fh.idx()
        valDT = aArea[nCell] / max(razon1)
        aDT.append(valDT)

    dt = CFL * min(aDT)                 # dt para la prox iteracion
    t += dt
    #print(t)
    aTime.append(t)

    hora += dt

    if hora >= 0.3:
        hora = 0


#####################################################################################################################################
# Exportar Datos
dfPM10 = pd.DataFrame(aUPM_timecell)
export = dfPM10.to_csv(r'visual/PM_estimadoTest.csv', index = None, header=True)

dfExac = pd.DataFrame(hExacta)
export = dfExac.to_csv(r'visual/PM_exactoTest.csv', index = None, header=True)

print('Listo!')
print(celdaPar)
