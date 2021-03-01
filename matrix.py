import numpy as np

matriz = np.matrix([[10,2,3], [3,4,5], [65,7,81]])

dett = np.linalg.det(matriz)


traspose = matriz.T

termino = traspose.item(-1)
apuntar = traspose.item((2,-1))


#OTRA MANERA DE ESCRIBIR UNA MATRIZ
lista = [[10,2,3], [3,4,5], [65,7,81]]

matriz = np.matrix(lista)

print(lista)
print(matriz)

# la matriz la converti en una Lista
nuevaLista = matriz.tolist()

print('Nueva Lista:', nuevaLista)

#creo la nueva fila a agregar
tupla = [1, 2, 10]
# la agrego al final de la matriz
nuevaLista.append(tupla)

print('Nueva Lista2:', nuevaLista)

nMatriz = np.matrix(nuevaLista)

print('Nueva Matriz con una nueva fila')
print(nMatriz)


#########################################################################################################################

## AGREGAR UNA NUEVA COLUMNA ##

# Crear una matriz con np.array o np.matrix es lo mismo
x = np.array([[10,20,30], [40,50,60]])

y = np.array([[100, 2], [200, 3]])

nMatriz2 = np.append(x, y, axis = 1)
print('Resultado')
print(nMatriz2)
# print(np.append(x, y, axis=1))

## AGREGAR UNA NUEVA FILA, USANDO VSTACK

nuevasFilas = [[1,2,3,4,5], [6,7,8,9,0]]
# Si las filas a agregar son mas que 1, entonces convertir a matris
nuevasFilas = np.array(nuevasFilas)
print(nuevasFilas)

matrizFinal = np.vstack((nMatriz2, nuevasFilas))

print(matrizFinal)