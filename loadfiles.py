import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
import numpy as np

##########################################################################################################
## PANDAS

# cargar un excel

df = pd.read_excel('centroides_alturas.xlsx', sheet_name=2, header = 0)

result = df.values.tolist()

for a in result:
    alturas = a


#cargar un CSV
data = pd.read_csv("PM_exacto.csv", sep = ",", encoding = "ISO=8859-1", decimal = ".")
# cargar un excel

df2 = pd.read_excel('asignaciones.xlsx', sheet_name=1, header = 0)
assig = df2.values.tolist()

for val in assig:
    asignation = val

def appoint():
    return asignation

def heights():
    return alturas