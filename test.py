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

arrayt = [2.54,3.12,4.445,5,6,2,3,1]
valorsillo = sum(arrayt)
arrayt = np.array(arrayt) * pow(10,-7)

fMovil_10 = [11.06, 2,3,4]


def moreless(val):
    if val == 0:
        return('It is equal to zaro')
    if val > 0:
        return('It is more than zero')
    if val < 0:
        return('It is less than zero')

valor = moreless(5)
print(valor)