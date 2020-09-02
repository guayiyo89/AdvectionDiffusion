<<<<<<< HEAD
import math

def valvorviento(vel):
    #recibo la velc del viento en m/s -> lo convierto a m/h
    velc = vel * 3600
    return vel

def calc_Vx(velc,dirv):
    # El dato de direccion recibido es de donde vino el viento.

    velc = valvorviento(velc)
    angulo = dirv

    if angulo >= 360:
        angulo = angulo - 360

    # angulo en radianes
    angulo = angulo * math.pi / 180

    vx = velc * math.sin(angulo)
    return vx

def calc_Vy(velc,dirv):
    # El dato de direccion recibido es de donde vino el viento.
    velc = valvorviento(velc)
    angulo = dirv

    if angulo >= 360:
        angulo = angulo - 360

    # angulo en radianes
    angulo = angulo * math.pi / 180

    vy = velc * math.cos(angulo)
=======
import math

def valvorviento(vel):
    #recibo la velc del viento en m/s -> lo convierto a m/h
    velc = vel * 3600
    return vel

def calc_Vx(velc,dirv):
    # El dato de direccion recibido es de donde vino el viento.

    velc = valvorviento(velc)
    angulo = dirv

    if angulo >= 360:
        angulo = angulo - 360

    # angulo en radianes
    angulo = angulo * math.pi / 180

    vx = velc * math.sin(angulo)
    return vx

def calc_Vy(velc,dirv):
    # El dato de direccion recibido es de donde vino el viento.
    velc = valvorviento(velc)
    angulo = dirv

    if angulo >= 360:
        angulo = angulo - 360

    # angulo en radianes
    angulo = angulo * math.pi / 180

    vy = velc * math.cos(angulo)
>>>>>>> 82192ae73ea8efbf3b6b90fd7171f97559717865
    return vy