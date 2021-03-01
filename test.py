import numpy as np

# USO DE NUMPY
valor = np.cos(np.pi)

print('El valor de cos de pi:', valor)


def miFuncion(num):
    doble = num * 2
    return doble

def multiplicar(num1, num2):
    resultado = num1 * num2
    return resultado

# DEFINIR UN ARRAY
arreglo = [] # esta vacio!

# UN CICLO WHILE
i = 1
while(i < 5):
    value = miFuncion(i)
    i = i + 1
    arreglo.append(value) # relleno el arreglo, con cada nuevo elemento al final de el.
    #print(arreglo)

#print('FINAL:', arreglo)
#print('El largo del arreglo es: ', len(arreglo))

""" fuunction Mifuncion(num){
    doble = num * 2;
} """ # FUNCION EN MATLAB

#def multiplicar(num1, num2)

# UN CICLO FOR

# FOR para recorrer un arreglo
for elemento in arreglo:
    result = multiplicar(elemento, 3)
    #print('EL valor es:', result)

# FOR UTILIZANDO UN RANGO
# for i in range(10):
    # print(i)

# range(partida, final, salto)
for valor in range(1, 11, 3):
    result = valor * 0.5
    # print(result)

# Definir un arreglo usando linspace
aValues = np.linspace(0, 10, 101)


print('aValues:', aValues)

primero = aValues[0]
ultimo = aValues[-2]
print('el primero es:', primero)
print('el ultimo es:', ultimo)