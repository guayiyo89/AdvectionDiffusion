import math

def calc_Vx(velc,dirv):
    # El dato de direccion recibido es de donde vino el viento.
    # Utilizaremos hacia donde se dirige.
    angulo = dirv + 180

    if angulo >= 360:
        angulo = angulo - 360

    # angulo en radianes
    angulo = angulo * math.pi / 180

    vx = velc * math.sin(angulo)
    return vx

def calc_Vy(velc,dirv):
    # El dato de direccion recibido es de donde vino el viento.
    # Utilizaremos hacia donde se dirige.
    angulo = dirv + 180

    if angulo >= 360:
        angulo = angulo - 360

    # angulo en radianes
    angulo = angulo * math.pi / 180

    vy = velc * math.cos(angulo)
    return vy
