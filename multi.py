import numpy as np
from multiprocessing import Process
import os
import winsound


#defino una funcion
def funcion(numero):
    print(os.getpid())
    for n in range(10):
        valor = n*n + n
        print(valor, "----->", numero)
        winsound.Beep(1500,120)

#para la ejecucion del codigo
if __name__ == '__main__':
    # lista para los procesos
    procesos = []
    # obtenemos la cantidad de nucleos
    nCores = os.cpu_count()
    print('El numero de nucleos es igual a: ', nCores)
    # generamos las instancias para los procesos
    for n in range(nCores):
        proceso = Process(target=funcion, args=(n,))
        # se anade a la lista de procesos
        procesos.append(proceso)

    print('Ejecutar...')

    for proceso in procesos:
        proceso.start()

    print('Espera...')
    for proceso in procesos:
        proceso.join()

    print('Vuelta al inicio...')