import random as rnd


def direction(value):
    if value >= 0 and value <= 45:
        coeff = rnd.random()
        if coeff >= 0.0 and coeff < 0.27:
            value = rnd.uniform(0,45)
        if coeff >= 0.27 and coeff < 0.44:
            value = rnd.uniform(45,90)
        if coeff >= 0.44 and coeff < 0.48:
            value = rnd.uniform(90,135)
        if coeff >= 0.48 and coeff < 0.5:
            value = rnd.uniform(135,180)
        if coeff >= 0.5 and coeff < 0.64:
            value = rnd.uniform(180,225)
        if coeff >= 0.64 and coeff < 0.86:
            value = rnd.uniform(225,270)
        if coeff >= 0.86 and coeff < 0.91:
            value = rnd.uniform(270,315)
        if coeff >= 0.91 and coeff <= 1.0:
            value = rnd.uniform(315,360)

    if value > 45 and value <= 90:
        coeff = rnd.random()
        if coeff >= 0.0 and coeff < 0.23:
            value = rnd.uniform(0,45)
        if coeff >= 0.23 and coeff < 0.45:
            value = rnd.uniform(45,90)
        if coeff >= 0.45 and coeff < 0.53:
            value = rnd.uniform(90,135)
        if coeff >= 0.53 and coeff < 0.56:
            value = rnd.uniform(135,180)
        if coeff >= 0.56 and coeff < 0.72:
            value = rnd.uniform(180,225)
        if coeff >= 0.72 and coeff < 0.93:
            value = rnd.uniform(225,270)
        if coeff >= 0.93 and coeff < 0.97:
            value = rnd.uniform(270,315)
        if coeff >= 0.97 and coeff <= 1.0:
            value = rnd.uniform(315,360)

    if value > 90 and value <= 135:
        coeff = rnd.random()
        if coeff >= 0.0 and coeff < 0.16:
            value = rnd.uniform(0,45)
        if coeff >= 0.16 and coeff < 0.28:
            value = rnd.uniform(45,90)
        if coeff >= 0.28 and coeff < 0.39:
            value = rnd.uniform(90,135)
        if coeff >= 0.39 and coeff < 0.48:
            value = rnd.uniform(135,180)
        if coeff >= 0.48 and coeff < 0.71:
            value = rnd.uniform(180,225)
        if coeff >= 0.71 and coeff < 0.94:
            value = rnd.uniform(225,270)
        if coeff >= 0.94 and coeff < 0.99:
            value = rnd.uniform(270,315)
        if coeff >= 0.99 and coeff <= 1.0:
            value = rnd.uniform(315,360)

    if value > 135 and value <= 180:
        coeff = rnd.random()
        if coeff >= 0.0 and coeff < 0.13:
            value = rnd.uniform(0,45)
        if coeff >= 0.13 and coeff < 0.18:
            value = rnd.uniform(45,90)
        if coeff >= 0.18 and coeff < 0.26:
            value = rnd.uniform(90,135)
        if coeff >= 0.26 and coeff < 0.39:
            value = rnd.uniform(135,180)
        if coeff >= 0.39 and coeff < 0.65:
            value = rnd.uniform(180,225)
        if coeff >= 0.65 and coeff < 0.90:
            value = rnd.uniform(225,270)
        if coeff >= 0.90 and coeff < 0.98:
            value = rnd.uniform(270,315)
        if coeff >= 0.98 and coeff <= 1.0:
            value = rnd.uniform(315,360)

    if value > 180 and value <= 225:
        coeff = rnd.random()
        if coeff >= 0.0 and coeff < 0.13:
            value = rnd.uniform(0,45)
        if coeff >= 0.13 and coeff < 0.18:
            value = rnd.uniform(45,90)
        if coeff >= 0.18 and coeff < 0.21:
            value = rnd.uniform(90,135)
        if coeff >= 0.21 and coeff < 0.27:
            value = rnd.uniform(135,180)
        if coeff >= 0.27 and coeff < 0.55:
            value = rnd.uniform(180,225)
        if coeff >= 0.55 and coeff < 0.83:
            value = rnd.uniform(225,270)
        if coeff >= 0.83 and coeff < 0.97:
            value = rnd.uniform(270,315)
        if coeff >= 0.97 and coeff <= 1.0:
            value = rnd.uniform(315,360)

    if value > 225 and value <= 270:
        coeff = rnd.random()
        if coeff >= 0.0 and coeff < 0.12:
            value = rnd.uniform(0,45)
        if coeff >= 0.12 and coeff < 0.18:
            value = rnd.uniform(45,90)
        if coeff >= 0.18 and coeff < 0.20:
            value = rnd.uniform(90,135)
        if coeff >= 0.20 and coeff < 0.22:
            value = rnd.uniform(135,180)
        if coeff >= 0.22 and coeff < 0.47:
            value = rnd.uniform(180,225)
        if coeff >= 0.47 and coeff < 0.79:
            value = rnd.uniform(225,270)
        if coeff >= 0.79 and coeff < 0.95:
            value = rnd.uniform(270,315)
        if coeff >= 0.95 and coeff <= 1.0:
            value = rnd.uniform(315,360)

    if value > 270 and value <= 315:
        coeff = rnd.random()
        if coeff >= 0.0 and coeff < 0.16:
            value = rnd.uniform(0,45)
        if coeff >= 0.16 and coeff < 0.22:
            value = rnd.uniform(45,90)
        if coeff >= 0.22 and coeff < 0.23:
            value = rnd.uniform(90,135)
        if coeff >= 0.23 and coeff < 0.24:
            value = rnd.uniform(135,180)
        if coeff >= 0.24 and coeff < 0.42:
            value = rnd.uniform(180,225)
        if coeff >= 0.42 and coeff < 0.70:
            value = rnd.uniform(225,270)
        if coeff >= 0.70 and coeff < 0.93:
            value = rnd.uniform(270,315)
        if coeff >= 0.93 and coeff <= 1.0:
            value = rnd.uniform(315,360)

    if value > 315 and value <= 360:
        coeff = rnd.random()
        if coeff >= 0.0 and coeff < 0.26:
            value = rnd.uniform(0,45)
        if coeff >= 0.26 and coeff < 0.34:
            value = rnd.uniform(45,90)
        if coeff >= 0.34 and coeff < 0.35:
            value = rnd.uniform(90,135)
        if coeff >= 0.35 and coeff < 0.36:
            value = rnd.uniform(135,180)
        if coeff >= 0.36 and coeff < 0.50:
            value = rnd.uniform(180,225)
        if coeff >= 0.50 and coeff < 0.73:
            value = rnd.uniform(225,270)
        if coeff >= 0.73 and coeff < 0.86:
            value = rnd.uniform(270,315)
        if coeff >= 0.86 and coeff <= 1.0:
            value = rnd.uniform(315,360)
    
    return value