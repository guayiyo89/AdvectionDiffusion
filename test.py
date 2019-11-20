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
print(sum(arrayt))