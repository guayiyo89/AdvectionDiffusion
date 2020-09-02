import openmesh as om
import numpy as np

mesh = om.TriMesh()

vh0 = mesh.add_vertex([-2, 4, 0])
vh1 = mesh.add_vertex([0, 2, 0])
vh2 = mesh.add_vertex([3, 4, 0])
vh3 = mesh.add_vertex([-2, 0, 0])
vh4 = mesh.add_vertex([3, 0, 0])

# esta seria la malla formada por los puntos.
#  0 ==== 2
#  |\  0 /|
#  | \  / |
#  |2  1 3|
#  | /  \ |
#  |/  1 \|
#  3 ==== 4

fh0 = mesh.add_face(vh0, vh1, vh2)
fh1 = mesh.add_face(vh1, vh3, vh4)
fh2 = mesh.add_face(vh0, vh3, vh1)

vh_list = [vh2, vh1, vh4]
fh3 = mesh.add_face(vh_list)

# write and read meshes
om.write_mesh('test.off', mesh)

# read the mesh
# mesh_2 = om.read_trimesh('test.off')

#inicializo arreglo para guardar longitud de lados
aLongEdge = np.zeros(8)

#arreglo que guarda valor de las areas
aArea = []

#arreglo q guarda los baricentros de celda
aCentCellx = []
aCentCelly = []

#arreglo que guarda los centros de interfaces
acenter_cell_x = np.zeros(8)
acenter_cell_y = np.zeros(8)

nCaras = 0



    for vh in mesh.fv(fh):
        indexv = vh.idx()
        punto = mesh.point(vh)
        vX = punto[0]
        vY = punto[1]
        verX.append(vX)
        verY.append(vY)


    # Calcularemos la longitud de cada lado de la celda
    i = 0
    for eh in mesh.fe(fh):
        nEdge = eh.idx()
        ax = verX[i - 1]
        bx = verX[i]
        ay = verY[i - 1]
        by = verY[i]
        valorX = pow(ax - bx, 2)
        valorY = pow(ay - by, 2)
        longitud = np.sqrt(valorX + valorY)

        #almaceno la longitud en un array
        aLongEdge[nEdge] = longitud
        # centro de cada lado de celda
        ncenter_cell_x = (ax + bx) / 2
        ncenter_cell_y = (ay + by) / 2

        acenter_cell_x[nEdge] = ncenter_cell_x
        acenter_cell_y[nEdge] = ncenter_cell_y

        i = i + 1

    # Calculen centro de las celdas
    # Calcular area de las celdas
    area = (verX[0] * (verY[1] - verY[2]) + verX[1] * (verY[2] - verY[0]) + verX[2] * (verY[0] - verY[1])) / 2
    aArea.append(np.abs(area))

    #calculo el baricentro de la celda
    nCentCellx = sum(verX) / 3
    nCentCelly = sum(verY) / 3

    aCentCellx.append(nCentCellx)
    aCentCelly.append(nCentCelly)

    nCaras = nCaras + 1

#############################################################################################################
#############################################################################################################

link_cell_edge = np.zeros((nCaras,3))
link_cell_cell = np.zeros((nCaras,3))

# obtengo los maximos y minimos de valores de los centros de interface
# asi reconozco que son frontera (solo para mallas cuyo limite es un rectangulo)
nMax_cx = np.max(acenter_cell_x)
nMax_cy = np.max(acenter_cell_y)
nMin_cx = np.min(acenter_cell_x)
nMin_cy = np.min(acenter_cell_y)


# Link Cell to Edge and Cell to Cell
for fh in mesh.faces():
    nFace = fh.idx()
    array_flag = []
    index = 0
    for eh in mesh.fe(fh):
        nEdge = eh.idx()
        link_cell_edge[nFace,index] = nEdge
        index = index + 1

        if acenter_cell_x[nEdge] == nMax_cx or acenter_cell_y[nEdge] == nMax_cy or acenter_cell_x[nEdge] == nMin_cx or acenter_cell_y[nEdge] == nMin_cy:
            array_flag.append(1)
            # Is a Frontier edge
        else:
            array_flag.append(0)
            # Doesn't is a frontier edge
    print(array_flag)

    # Variable to save the indexes of every neighbour face
    auxCell1 = []
    for fh in mesh.ff(fh):
        nFace1 = fh.idx()
        auxCell1.append(nFace1)

    index = 0
    aux_idx = 0 # this one works as a counter for auxCell1

    # This for loop goes trough the values of the array_flag
    for value in array_flag:
        if value == 1:
            link_cell_cell[nFace, index] = int(-1)
        else:
            link_cell_cell[nFace, index] = auxCell1[aux_idx]
            # this is for assure auxCell1 only will be goes trough when the edge
            aux_idx = aux_idx + 1
        index = index + 1




print('Feces: ',link_cell_cell)
print('Edges: ',link_cell_edge)


# esta seria la malla formada por los puntos.
#  0 ==== 2
#  |\  0 /|
#  | \  / |
#  |2  1 3|
#  | /  \ |
#  |/  1 \|
#  3 ==== 4
