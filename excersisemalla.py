<<<<<<< HEAD
import openmesh as om
import numpy as np

mesh = om.TriMesh()

vh0 = mesh.add_vertex([0, 23, 0])
vh1 = mesh.add_vertex([12, 23, 0])
vh2 = mesh.add_vertex([23, 23, 0])
vh3 = mesh.add_vertex([38, 23, 0])
vh4 = mesh.add_vertex([47, 23, 0])
vh5 = mesh.add_vertex([47, 17, 0])
vh6 = mesh.add_vertex([47, 7, 0])
vh7 = mesh.add_vertex([47, 0, 0])
vh8 = mesh.add_vertex([38, 0, 0])
vh9 = mesh.add_vertex([26, 0, 0])
vh10 = mesh.add_vertex([13, 0, 0])
vh11 = mesh.add_vertex([0, 0, 0])
vh12 = mesh.add_vertex([0, 5, 0])
vh13 = mesh.add_vertex([0, 16, 0])

vh14 = mesh.add_vertex([9, 18, 0])
vh15 = mesh.add_vertex([29, 17, 0])
vh16 = mesh.add_vertex([14, 11, 0])
vh17 = mesh.add_vertex([32, 8, 0])

# esta seria la malla formada por los puntos.
#  0 ==== 2
#  |\  0 /|
#  | \  / |
#  |2  1 3|
#  | /  \ |
#  |/  1 \|
#  3 ==== 4

fh0 = mesh.add_face(vh0, vh13, vh14)
fh1 = mesh.add_face(vh14, vh1, vh0)
fh2 = mesh.add_face(vh13, vh16, vh14)
fh3 = mesh.add_face(vh14, vh2, vh1)
fh4 = mesh.add_face(vh16, vh15, vh14)
fh5 = mesh.add_face(vh13, vh12, vh16)
fh6 = mesh.add_face(vh14, vh15, vh2)
fh7 = mesh.add_face(vh16, vh17, vh15)
fh8 = mesh.add_face(vh16, vh12, vh10)
fh9 = mesh.add_face(vh3, vh2, vh15)
fh10 = mesh.add_face(vh9, vh17, vh16)
fh11 = mesh.add_face(vh5, vh15, vh17)
fh12 = mesh.add_face(vh16, vh10, vh9)
fh13 = mesh.add_face(vh10, vh12, vh11)

fh14 = mesh.add_face(vh3, vh15, vh5)
fh15 = mesh.add_face(vh8, vh17, vh9)
fh16 = mesh.add_face(vh17, vh6, vh5)
fh17 = mesh.add_face(vh4, vh3, vh5)
fh18 = mesh.add_face(vh8, vh6, vh17)
fh19 = mesh.add_face(vh6, vh8, vh7)
# write and read meshes
om.write_mesh('test.off', mesh)



aFace_n = []
aEdge = []
for eh in mesh.fe(fh1):
    nEdge = eh.idx()
    aEdge.append(nEdge)
for fh in mesh.ff(fh1):
    nFace_n = fh.idx()
    aFace_n.append(nFace_n)

print('Edges: ',aEdge)
=======
import openmesh as om
import numpy as np

mesh = om.TriMesh()

vh0 = mesh.add_vertex([0, 23, 0])
vh1 = mesh.add_vertex([12, 23, 0])
vh2 = mesh.add_vertex([23, 23, 0])
vh3 = mesh.add_vertex([38, 23, 0])
vh4 = mesh.add_vertex([47, 23, 0])
vh5 = mesh.add_vertex([47, 17, 0])
vh6 = mesh.add_vertex([47, 7, 0])
vh7 = mesh.add_vertex([47, 0, 0])
vh8 = mesh.add_vertex([38, 0, 0])
vh9 = mesh.add_vertex([26, 0, 0])
vh10 = mesh.add_vertex([13, 0, 0])
vh11 = mesh.add_vertex([0, 0, 0])
vh12 = mesh.add_vertex([0, 5, 0])
vh13 = mesh.add_vertex([0, 16, 0])

vh14 = mesh.add_vertex([9, 18, 0])
vh15 = mesh.add_vertex([29, 17, 0])
vh16 = mesh.add_vertex([14, 11, 0])
vh17 = mesh.add_vertex([32, 8, 0])

# esta seria la malla formada por los puntos.
#  0 ==== 2
#  |\  0 /|
#  | \  / |
#  |2  1 3|
#  | /  \ |
#  |/  1 \|
#  3 ==== 4

fh0 = mesh.add_face(vh0, vh13, vh14)
fh1 = mesh.add_face(vh14, vh1, vh0)
fh2 = mesh.add_face(vh13, vh16, vh14)
fh3 = mesh.add_face(vh14, vh2, vh1)
fh4 = mesh.add_face(vh16, vh15, vh14)
fh5 = mesh.add_face(vh13, vh12, vh16)
fh6 = mesh.add_face(vh14, vh15, vh2)
fh7 = mesh.add_face(vh16, vh17, vh15)
fh8 = mesh.add_face(vh16, vh12, vh10)
fh9 = mesh.add_face(vh3, vh2, vh15)
fh10 = mesh.add_face(vh9, vh17, vh16)
fh11 = mesh.add_face(vh5, vh15, vh17)
fh12 = mesh.add_face(vh16, vh10, vh9)
fh13 = mesh.add_face(vh10, vh12, vh11)

fh14 = mesh.add_face(vh3, vh15, vh5)
fh15 = mesh.add_face(vh8, vh17, vh9)
fh16 = mesh.add_face(vh17, vh6, vh5)
fh17 = mesh.add_face(vh4, vh3, vh5)
fh18 = mesh.add_face(vh8, vh6, vh17)
fh19 = mesh.add_face(vh6, vh8, vh7)
# write and read meshes
om.write_mesh('test.off', mesh)



aFace_n = []
aEdge = []
for eh in mesh.fe(fh1):
    nEdge = eh.idx()
    aEdge.append(nEdge)
for fh in mesh.ff(fh1):
    nFace_n = fh.idx()
    aFace_n.append(nFace_n)

print('Edges: ',aEdge)
>>>>>>> 82192ae73ea8efbf3b6b90fd7171f97559717865
print('Faces: ',aFace_n)