import numpy as np
import matplotlib.pyplot as plt

def fun(x,y):
    #aqui colocamos f(x,y)
    nf = 4 * np.exp(0.8 * x) - 0.5 * y
    return nf

aY = [2]                     #Defino la condicion inicial y para x = 0

aX = np.linspace(0,4,9)     #Defino los valores de x
h = aX[1] - aX[0]           #calculo el tamano del paso


for i in range(len(aX)-1):
    x = aX[i]
    y = aY[i]

    #evaluo en dy/dx
    k1 = fun(aX[i], aY[i-1])
    k2 = fun(aX[i] + 0.5 * h, aY[i-1] + 0.5 * h * k1)
    k3 = fun(aX[i] + 0.5 * h, aY[i-1] + 0.5 * h * k2)
    k4 = fun(aX[i] + h, aY[i-1] + h * k3)

    print(k1, k2, k3, k4)

    y = aY[i-1] + (1/6) * h * (k1 + 2*k2 + 2*k3 + k4)

    aY.append(y)

print(aX, aY)

### GRAFICA
plt.plot(aX,aY)
plt.show()