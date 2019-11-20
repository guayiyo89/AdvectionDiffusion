import pandas as pd
import numpy as np

# Velc viento
data = pd.read_csv(".\datos\datos_velviento_encinas.csv", sep = ";", encoding = "ISO=8859-1", decimal = ",")
result = data.values.tolist()
velv = []

for a in result:
    if pd.notna(a[2]):
        valor = a[2] * 3.6
        velv.append(valor)

valores1 = []; valores2 = []; valores3 = []; valores4 = []; valores5 = []; valores6 = []; valores7 = []
for value in velv:
    if value >= 0 and value < 3:
        valores1.append(value)
    if value >= 3 and value < 6:
        valores2.append(value)
    if value >= 6 and value < 9:
        valores3.append(value)
    if value >= 9 and value < 12:
        valores4.append(value)
    if value >= 12 and value < 15:
        valores5.append(value)
    if value >= 15 and value < 20:
        valores6.append(value)
    if value > 20:
        valores7.append(value)

print(len(valores1)/len(velv)); print(len(valores2)/len(velv)); print(len(valores3)/len(velv))
print(len(valores4)/len(velv)); print(len(valores5)/len(velv)); print(len(valores6)/len(velv)); print(len(valores7)/len(velv))