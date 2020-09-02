<<<<<<< HEAD
import numpy as np
import pandas as pd

# Direction viento
data = pd.read_csv("ambientales_var\PM10.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
result = data.values.tolist()
PM10 = []

for a in result:
    if pd.notna(a[2]):
        PM10.append(a[2])

data = pd.read_csv("ambientales_var\PM25.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
result = data.values.tolist()
PM25 = []

for a in result:
    if pd.notna(a[2]):
        PM25.append(a[2])

print('LARGOS')
print(len(PM10))
print(len(PM25))

sPM10 = []
sPM25 = []

for n in PM10:
    if n > 50:
        sPM10.append(n)

for n in PM25:
    if n > 25:
        sPM25.append(n)

print('LARGOS')
print(len(sPM10))
=======
import numpy as np
import pandas as pd

# Direction viento
data = pd.read_csv("ambientales_var\PM10.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
result = data.values.tolist()
PM10 = []

for a in result:
    if pd.notna(a[2]):
        PM10.append(a[2])

data = pd.read_csv("ambientales_var\PM25.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
result = data.values.tolist()
PM25 = []

for a in result:
    if pd.notna(a[2]):
        PM25.append(a[2])

print('LARGOS')
print(len(PM10))
print(len(PM25))

sPM10 = []
sPM25 = []

for n in PM10:
    if n > 50:
        sPM10.append(n)

for n in PM25:
    if n > 25:
        sPM25.append(n)

print('LARGOS')
print(len(sPM10))
>>>>>>> 82192ae73ea8efbf3b6b90fd7171f97559717865
print(len(sPM25))