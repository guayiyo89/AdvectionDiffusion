import numpy as np
import openmesh as om
import random
n = 4
m = 3
#a = n*m
matriz = []

for i in range(n):
    matriz.append([])
    for j in range(m):
        matriz[i].append(random.randint(0, 100))

dime = type(-1.0)

arrayt = [2,3,4,5,6,2,3,10]
arrayt = np.array(arrayt) * pow(10,-7)
print(arrayt)

fMovil_10 = [11.06, 2,3,4]
print(fMovil_10)