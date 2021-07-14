import numpy as np
import openmesh as om
import pandas as pd

mesh = om.TriMesh()

# vertices de malla
vh0 = mesh.add_vertex([0, 0, 0])
vh1 = mesh.add_vertex([2, 2, 0])
vh2 = mesh.add_vertex([0, 4, 0])
vh3 = mesh.add_vertex([4, 0, 0])
vh4 = mesh.add_vertex([4, 4, 0])

fh0 = mesh.add_face(vh0, vh1, vh2)
fh1 = mesh.add_face(vh2, vh1, vh4)
fh2 = mesh.add_face(vh4, vh1, vh3)
fh3 = mesh.add_face(vh3, vh1, vh0)

om.write_mesh('example.off', mesh)