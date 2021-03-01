import numpy as np
import openmesh as om
import pandas as pd

#### Defino el tipo de Malla
mesh = om.TriMesh() # TriMesh es una funcion de OpenMesh (om)

#Al momento de crear los vertices se agrega un Arreglo  [x,y,z] (z=0 para asegurar 2D)
#a la funcion mesh.add_vertex()
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

#Generar archivo
om.write_mesh('test.off', mesh)

##### COMANDOS BASICOS

#Recorro todo los nodo de la malla
for vh in mesh.vertices():
    #IMPRIMIR LAS COORDENADAS DE CADA VERTICE
    punto = mesh.point(vh)
    # print(punto)

#Recorrer todas las interfaces
contador = 0
for eh in mesh.edges():
    eIdx = eh.idx() # devuelve el indice global del elemento
    contador = contador + 1

eIdx = eIdx + 1

print("total interfaces: ",contador) # valor del total de numero de interfaces de la malla

#arreglo q guardaba la longitud de los lados
aLongLados = np.zeros(contador)

#arreglos que guardan los centros de cada lado
aCenint_x  = np.zeros(contador)
aCenint_y  = np.zeros(contador)

#arreglo guarda los baricentros de la celda
aBarX = []
aBarY = []
#arreglo q guarda las aresa
aArea = []

#arreglos q guardasn las normales
aNormal_xf = []
aNormal_yf = []

#Distancia centro celda a centro lado
aDist_cell_edge = []

# sumatoria de los inversos de las distancias nodo - centro (vertices)
aSum_dist_node_cell = []

# dictancia nodo - centro de celda
aDist_node_cell = np.zeros((8,8))

# relacion celda a sus bordes
link_cell_to_edge = np.zeros((66,3))
link_edge_to_cell = []
# relacion celda a celda vecinas
link_cell_to_cell = np.zeros((66,3))
# relacion celda vertices
link_cell_to_vertex = np.zeros((66,3))

#RECORRER
for fh in mesh.faces():
    #IMPRIMIR LOS INDICES DE CADA CELDA(faces)
    index = fh.idx() # devuelve el indice de la celda
    #print("Cara: ",index)
    #Recorrer los vertices de cada cara y mostrar sus indices
    for vh in mesh.fv(fh):
        idx = vh.idx() # devuelve indice del VERTICE
        #print(idx)
    for fh in mesh.ff(fh):
        nCell = fh.idx()
        #print(nCell)
    for eh in mesh.fe(fh):
        nEdge = eh.idx()
        #print('Lados',nEdge)

#################################################################
#CALCULAR LONGITUD DE LOS LADOS DE LAS CELDAS


for fh in mesh.faces():
    nCell = fh.idx()

    aVetX = [] # guardar las coordenadas X de los vertices
    aVetY = [] # guardar las coordenadas Y de los vertices
    for vh in mesh.fv(fh):
        punto = mesh.point(vh) # devolver un arreglo de long 3
        vetX = punto[0]
        vetY = punto[1]
        aVetX.append(vetX)
        aVetY.append(vetY)
    # aVetX y aVetY deberian tener cada una de las coordenadas de los vertices

    # Calcular el baricentro de la celda
    barX = sum(aVetX)/3
    barY = sum(aVetY)/3

    aBarX.append(barX)
    aBarY.append(barY)

    # AREA DE LAS CELDAS
    nArea = aVetX[0] * (aVetY[2] - aVetY[1]) + aVetX[1] * (aVetY[0] - aVetY[2]) + aVetX[2] * (aVetY[1] - aVetY[0])
    nArea = 0.5 * np.absolute(nArea)
    aArea.append(nArea)

    #auxiliar para el calculo distancia centro celda - centro lado
    aux_dist_cCell_cEdge = []

    iaux = 0 # contador auxiliar

    #auxiliar q guarda las normales (largo 3)
    auxNormal_x = []
    auxNormal_y = []

    for eh in mesh.fe(fh):
        # CALCULAR LA DISTANCIA
        nEdge = eh.idx()
        distX = aVetX[iaux -1] - aVetX[iaux]
        distY = aVetY[iaux -1] - aVetY[iaux]

        # calculado el punto medio de cada lado
        nCentLadoX = 0.5 * (aVetX[iaux - 1] + aVetX[iaux])
        nCentLadoY = 0.5 * (aVetY[iaux - 1] + aVetY[iaux])
        #lo guardamos en el arreglo
        aCenint_x[nEdge] = nCentLadoX
        aCenint_y[nEdge] = nCentLadoY

        #formula de la distancia (Longitud segmento)
        distancia = np.sqrt(pow(distX,2) + pow(distY,2))
        aLongLados[nEdge] = distancia # guardamos en el arreglo

        # Calcular la distancia baricentro a centro de segmento
        dif_x = aBarX[nCell] - nCentLadoX
        dif_y = aBarY[nCell] - nCentLadoY
        distanceCL = np.sqrt(pow(dif_x,2) + pow(dif_y,2))

        aux_dist_cCell_cEdge.append(distanceCL)

        # tangente unitaria
        t_xf = -distX / distancia
        t_yf = -distY / distancia
        # normal unitaria
        auxNormal_x.append(t_yf)
        auxNormal_y.append(-t_xf)
        #Haganla!
        iaux = iaux + 1
    

    # GUARDO LAS 3 NORMALES DE CADA CARA
    aNormal_xf.append(auxNormal_x)
    aNormal_yf.append(auxNormal_y)

    #Las distancias centro a centro se guardan
    aDist_cell_edge.append(aux_dist_cCell_cEdge)

#obtenemos los maximos y los minimos de los centros de los bordes
nMax_Cx = np.amax(aCenint_x)
nMax_Cy = np.amax(aCenint_y)
nMin_Cx = np.amin(aCenint_x)
nMin_Cy = np.amin(aCenint_y)

#CALCULOS DE LA ECUACION
# tamanp de vector de tiempo
sizeT = 20
# Arreglo de tiempo
#aTime = np.linspace(0,sizeT,np.size(temperatura))

dt = 1 / sizeT

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
        #print("celda: ",nCell,fronteraux)
        index = index + 1
    print("celda: ",nCell,fronteraux)
    

    for fh in mesh.ff(fh):
        nCell_v = fh.idx()
        auxCell.append(nCell_v)
    
    # OBTENEMOS EL VECTOR DE LINK DE LA CELDA CON SUS VECINOS
    index = 0
    auxidx = 0 # recorre el auxiliar auxCell
    for value in fronteraux:
        if value == 1:
            aux = int(-1)
            link_cell_to_cell[nCell, index] = aux
        else:
            link_cell_to_cell[nCell, index] = auxCell[auxidx]
            auxidx = auxidx + 1
        #print(link_cell_to_cell)
        index = index + 1

    # OBTENEMOS LA MATRIZ DE LINK DE LA CELDA CON SUS VERTICES
    index = 0
    for vh in mesh.fv(fh):
        nVertex = vh.idx()
        link_cell_to_vertex[nCell,index] = nVertex
        index = index + 1

aW_f = []

#######################################################################################
#######################################################################################
#DEFINO LA CONCEBTRACION INICIAL

historico = []
historico.append(aU)
#Defino concentracion en el borde
aF = 0 #cte
#Difusion
Diff = 1
#######################################################################################
#######################################################################################

for t in range(sizeT):
    # auxiliar para el calculo de la concentracion
    auxU = []

    for fh in mesh.faces():
        # obtenemos indice de la celda
        nCell = fh.idx()

        # Almacenamos los flujos de la cara
        aFlux = []

        #recorrer los lados
        for i in range(3):
            # obtengo el indice del lado y de la celda vecina
            nEdge = int(link_cell_to_edge[nCell, i])
            nCell_v = int(link_cell_to_cell[nCell, i])
            # obtengo los indices de los vertices
            a = aVetX[i]
            b = aVetY[i-1]

            if nCell_v == -1:
                #Calculo la distancia entre el baricentro y el centro del borde de frontera
                distx = aBarX[nCell] - aCenint_x[nEdge]
                disty = aBarY[nCell] - aCenint_y[nEdge]
                l_f = np.sqrt(pow(distx,2) + pow(disty,2))
                # Calculo la distancia delta como productro de la dist entre los baricentros por la normal del segmento
                delta_distx = aNormal_xf[nCell][i] * distx
                delta_disty = aNormal_yf[nCell][i] * disty
                delta_dist = np.sqrt(pow(delta_distx,2) + pow(delta_disty,2))
                # Calculo el flujo
                flux = ((aF - aU[nCell]) / delta_dist) * aLongLados[nEdge] * Diff
                aFlux.append(flux)

            else:
                # Calculo la distancia entre los baricentros
                distx = aBarX[nCell] - aBarX[nCell_v]
                disty = aBarY[nCell] - aBarY[nCell_v]
                l_f = np.sqrt(pow(distx,2) + pow(disty,2))
                # Calculo la distancia delta como productro de la dist entre los baricentros por la normal del segmento
                delta_distx = aNormal_xf[nCell][i] * distx
                delta_disty = aNormal_yf[nCell][i] * disty
                delta_dist = np.sqrt(pow(delta_distx,2) + pow(delta_disty,2))
                # Calculo el flujo
                flux = ((aU[nCell_v] - aU[nCell]) / delta_dist) * aLongLados[nEdge] * Diff
                aFlux.append(flux)

        Ucell = aU[nCell] + (dt * sum(aFlux)) / aArea[nCell] # calculo concentracion en la celda

        # adjunto el vector con los 8 valores de U para las celdas
        auxU.append(Ucell)
        #Fin fh
        
        
    # agrego al historico
    aU = auxU
    historico.append(aU)


#Creo el CSV

dfU = pd.DataFrame(historico)
#print(dfU)

#exportamos un CSV
export = dfU.to_csv(r'datosU.csv', index = None, header=True)
