import numpy as np
import pandas as pd

# Direction viento
data = pd.read_csv("\PM10.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
result = data.values.tolist()
PM10 = []

for a in result:
    if pd.notna(a[2]):
        PM10.append(a[2])

data = pd.read_csv("PM25.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
result = data.values.tolist()
PM25 = []

for a in result:
    if pd.notna(a[2]):
        PM25.append(a[2])

print(PM10)