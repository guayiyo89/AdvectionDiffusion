import numpy as np
import openmesh as om
import pandas as pd
import exactas as fnc
import random

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

om.write_mesh('testmesh.off', mesh)