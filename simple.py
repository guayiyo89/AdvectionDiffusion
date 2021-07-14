import numpy as np
import openmesh as om
import pandas as pd
import exactas as fnc
import random

#### Defino el tipo de Malla
mesh = om.TriMesh() # TriMesh es una funcion de OpenMesh (om)


# vertices de malla
vh0 = mesh.add_vertex([-1, 1, 0])
vh1 = mesh.add_vertex([1, 1, 0])
vh2 = mesh.add_vertex([-1, -1, 0])
vh3 = mesh.add_vertex([1, -1, 0])

fh0 = mesh.add_face(vh1, vh0, vh3)
fh1 = mesh.add_face(vh3, vh0, vh2)

om.write_mesh('simple.off', mesh)