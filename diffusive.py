import math
import random as rd


# Boltzmann
kB = 1.38064 * pow(10, -23)
# Presion (Pa)
# los valores de prsion varian aleatoriamente entre estos valores
P = rd.uniform(100500, 102700)
# Viscosidad del aire  por temperatura (tabla)
nu_10 = 1.67 * pow(10,-5)
nu00 = 1.72 * pow(10,-5)
nu10 = 1.76 * pow(10,-5)
nu20 = 1.81 * pow(10,-5)
nu30 = 1.86 * pow(10,-5)
nu40 = 1.91 * pow(10,-5)

# Diametro moleculas de aire
dm = 4 * pow(10,-10)

#VISCOSIDAD AIRE (cP (Pa s2))
def visc_din(temp):
    if temp <= 0.0:
        val1 = -10; val2 = 0.0
        fval1 = nu_10; fval2 = nu00
    if temp >= 0.0 and temp < 10.0:
        val1 = 0; val2 = 10
        fval1 = nu00; fval2 = nu10
    if temp >= 10.0 and temp < 20.0:
        val1 = 10.0; val2 = 20.0
        fval1 = nu10; fval2 = nu20
    if temp >= 20.0 and temp < 30.0:
        val1 = 0.0; val2 = 10.0
        fval1 = nu20; fval2 = nu30
    if temp >= 30.0 and temp < 40.0:
        val1 = 30.0; val2 = 40.0
        fval1 = nu30; fval2 = nu40

    nu = fval1 + ((fval2 - fval1)/(val2 - val1)) * (temp - val1)

    return nu


def coeffdif(temp, diam, lamb, nu):
    diam = diam * pow(10,-6) # recibo el diametro
    value = -1.10 * diam / lamb
    cC = 1 + lamb / diam * (1.257 + 0.4 * math.exp(value))
    nu = visc_din(temp)
    temp = temp + 273.15
    diff = (2 * kB * temp) / (3 * math.pi * nu * diam / 2)
    return diff

def freepath(temp):
    temp = temp + 273.15
    lamb = (kB * temp) / (4 * math.sqrt(2) * math.pi * pow(dm,2) * P)
    return lamb

def cunningham(diam,temp):
    temp = temp + 273.15
    lamb = freepath(temp)
    value = -1.10 * diam / lamb
    cCoef = 1 + 2 * lamb / diam * (1.257 + 0.4 * math.exp(value))
    return cCoef

#CONVERTIR VIENTO A METRO POR SEGUNDO Y DE SALIDA RECONVERTIR
def vdep(temp, diam, uC, viento):
    diam = diam * pow(10,-6) # recibo el diametro
    cCoef = cunningham(diam,temp)
    viento = viento / 1000
    g = 9.81 #aceleracion m/s2
    # time particle relax
    tau = uC * cCoef * pow(diam,2) / (18 * visc_din(temp))
    # von Karmann Constant
    k = 0.4
    # Drag Coeff
    Cd = 1.3 * pow(10,-3)
    Vc1 = 1/(1-k) * Cd * viento

    #alpha (m3)
    ap = 100
    # Schmidt Number
    Sc = visc_din(temp)/diam

    viento = viento * 0.27778

    # Stoke number (dn = 1)
    St = tau * viento

    Vd1 = -ap * 0.1 + (1/k) * cCoef * viento * math.sqrt(Sc) + St/1000
    Vg = tau * g

    VC = Vc1 + Vg
    VD = Vd1 + Vg

    Vdep = 1/VC + 1/VD + (Vg / (VC * VD))
    Vdep = 1/Vdep
    return Vdep


def pond(temp):
    if temp <= 0:
        value = 0.02
    if temp > 0 and temp <= 5:
        value = 0.00875
    if temp > 5 and temp <= 10:
        value = 0.006
    if temp > 10 and temp <= 15:
        value = 0.00333
    if temp > 15 and temp <= 20:
        value = 0.00125
    if temp > 20 and temp <= 25:
        value = 0.00075
    if temp > 25:
        value = 0.0005
    return value
