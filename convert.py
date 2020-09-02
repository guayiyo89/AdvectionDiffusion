import numpy as np
import pandas as pd

f = lambda x : (x.replace(",","."))
# PM 2.5
data = pd.read_csv("datos\datos_140101_190930.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
hola = data.head()

df = pd.DataFrame(data, columns = ['FECHA (YYMMDD)', 'HORA (HHMM)', 'Registros validados', 'Registros preliminares', 'Registros no validados'])


result = df.values.tolist()

datanew = np.asarray(result)


pm = []
index = 0
for a in result:
    c = -1 # default
    if pd.notna(a[4]):
        c = a[4]
    if pd.notna(a[3]):
        c = a[3]
    if pd.notna(a[2]):
        c = a[2]

    b = [a[0], c]
    pm.append(b)

df2 = pd.DataFrame(pm, columns = ['FECHA (YYMMDD)', 'Registros'])

export = df2.to_csv(r'nielolDia_25.csv', index = None, header=True)

