import multiprocessing

multiprocessing.set_start_method('spawn', True)     # para q funcione el Pool
from multiprocessing.pool import Pool

def cubo(n):
    valor = n **3
    return valor

def sum_numeros(num):
    s = 0
    for i in range(num + 1):
        s += i * i
    return s

if __name__ == '__main__':
    print(Pool().map(cubo, [2,3,5]))

    numeros = range(2000)
    p = Pool()          #agrupa los procesos
    result = p.map(sum_numeros, numeros)
    print(result)

    p.close()
    p.join()
