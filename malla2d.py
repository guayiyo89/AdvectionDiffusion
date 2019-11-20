import numpy as np
import openmesh as om

mesh = om.TriMesh()

# vertices de malla
vh0 = mesh.add_vertex([0, 11.35646])
vh1 = mesh.add_vertex([6.479, 11.35646])
vh2 = mesh.add_vertex([12.254, 11.35646])
vh3 = mesh.add_vertex([17.732, 11.35646])
vh4 = mesh.add_vertex([17.732, 4.744824])
vh5 = mesh.add_vertex([17.732, 0])
vh6 = mesh.add_vertex([8.712, 0])
vh7 = mesh.add_vertex([0, 0])
vh8 = mesh.add_vertex([0, 4.900392])
vh9 = mesh.add_vertex([0, 7.467264])
# vertices interiores
vh10 = mesh.add_vertex([3.905, 9.022944])
vh11 = mesh.add_vertex([8.129, 8.567352])
vh12 = mesh.add_vertex([5.808, 5.678232])
vh13 = mesh.add_vertex([10.692, 6.311616])
vh14 = mesh.add_vertex([13.607, 8.711808])
vh15 = mesh.add_vertex([15.136, 6.900552])
vh16 = mesh.add_vertex([3.553, 2.83356])
vh17 = mesh.add_vertex([7.568, 2.922456])
vh18 = mesh.add_vertex([13.024, 3.411384])

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

#iteramos sobre los vertices de la cara
for fh in mesh.faces():
    aVtx = []
    aVty = []
    for vh in mesh.fv(fh):
        point = mesh.point(vh)
        index = vh.idx()
        aVtx.append(point[0])
        aVty.append(point[1])
    # calculamos el baricentro de la celda
    sumax = aVtx[0] + aVtx[1] + aVtx[2]
    nCenterx = 0.5 * sumax
    aCenterx.append(nCenterx)
    sumay = aVty[0] + aVty[1] + aVty[2]
    nCentery = 0.5 * sumay
    aCentery.append(nCentery)
    # calculamos el area de la celda
    print(fh.idx())
    nArea = aVtx[0] * (aVty[2] - aVty[1]) + aVtx[1] * (aVty[0] - aVty[2]) + aVtx[2] * (aVty[1] - aVty[0])
    nArea = 0.5 * np.absolute(nArea)
    aArea.append(nArea)
print("Area")
print(aArea)
