import random as rnd

def valor_velc(value):
    if value >= 0 and value <= 3:
        coeff = rnd.random()
        if coeff >= 0.0 and coeff < 0.47:
            value = rnd.uniform(0,3)
        if coeff >= 0.47 and coeff < 0.84:
            value = rnd.uniform(3,6)
        if coeff >= 0.84 and coeff < 0.97:
            value = rnd.uniform(6,9)
        if coeff >= 0.97:
            value = rnd.uniform(9,12)

    if value > 3 and value <= 6:
        coeff = rnd.random()
        if coeff >= 0.0 and coeff < 0.35:
            value = rnd.uniform(0,3)
        if coeff >= 0.35 and coeff < 0.75:
            value = rnd.uniform(3,6)
        if coeff >= 0.75 and coeff < 0.95:
            value = rnd.uniform(6,9)
        if coeff >= 0.95 and coeff < 0.99:
            value = rnd.uniform(9,12)
        if coeff >= 0.99:
            value = rnd.uniform(12,15)

    if value > 6 and value <= 9:
        coeff = rnd.random()
        if coeff >= 0.0 and coeff < 0.25:
            value = rnd.uniform(0,3)
        if coeff >= 0.25 and coeff < 0.56:
            value = rnd.uniform(3,6)
        if coeff >= 0.56 and coeff < 0.86:
            value = rnd.uniform(6,9)
        if coeff >= 0.86 and coeff < 0.96:
            value = rnd.uniform(9,12)
        if coeff >= 0.96 and coeff < 0.99:
            value = rnd.uniform(12,15)
        if coeff >= 0.99:
            value = rnd.uniform(15,18)

    if value > 9 and value <= 12:
        coeff = rnd.random()
        if coeff >= 0.0 and coeff < 0.12:
            value = rnd.uniform(1,3)
        if coeff >= 0.12 and coeff < 0.37:
            value = rnd.uniform(3,6)
        if coeff >= 0.37 and coeff < 0.66:
            value = rnd.uniform(6,9)
        if coeff >= 0.66 and coeff < 0.87:
            value = rnd.uniform(9,12)
        if coeff >= 0.87 and coeff < 0.97:
            value = rnd.uniform(12,15)
        if coeff >= 0.97 and coeff < 0.99:
            value = rnd.uniform(15,20)
        if coeff >= 0.99:
            value = rnd.uniform(20,25)

    if value > 12 and value <= 15:
        coeff = rnd.random()
        if coeff >= 0.0 and coeff < 0.09:
            value = rnd.uniform(0,3)
        if coeff >= 0.09 and coeff < 0.28:
            value = rnd.uniform(3,6)
        if coeff >= 0.28 and coeff < 0.48:
            value = rnd.uniform(6,9)
        if coeff >= 0.48 and coeff < 0.76:
            value = rnd.uniform(9,12)
        if coeff >= 0.76 and coeff < 0.94:
            value = rnd.uniform(12,15)
        if coeff >= 0.94 and coeff < 0.99:
            value = rnd.uniform(15,20)
        if coeff >= 0.99:
            value = rnd.uniform(20,30)

    if value > 15 and value <= 20:
        coeff = rnd.random()
        if coeff >= 0.0 and coeff < 0.01:
            value = rnd.uniform(0,3)
        if coeff >= 0.01 and coeff < 0.12:
            value = rnd.uniform(3,6)
        if coeff >= 0.12 and coeff < 0.27:
            value = rnd.uniform(6,9)
        if coeff >= 0.27 and coeff < 0.50:
            value = rnd.uniform(9,12)
        if coeff >= 0.50 and coeff < 0.79:
            value = rnd.uniform(12,15)
        if coeff >= 0.79 and coeff < 0.96:
            value = rnd.uniform(15,20)
        if coeff >= 0.96:
            value = rnd.uniform(20,40)

    if value > 20:
        coeff = rnd.random()
        if coeff >= 0.0 and coeff < 0.05:
            value = rnd.uniform(0,6)
        if coeff >= 0.05 and coeff < 0.15:
            value = rnd.uniform(6,9)
        if coeff >= 0.15 and coeff < 0.30:
            value = rnd.uniform(9,12)
        if coeff >= 0.30 and coeff < 0.63:
            value = rnd.uniform(12,15)
        if coeff >= 0.63 and coeff < 0.93:
            value = rnd.uniform(15,20)
        if coeff >= 0.93:
            value = rnd.uniform(20,60)

    return value

def temperatura(value):
    if value < 0:
        coeff = rnd.random()
        if coeff >= 0 and coeff < 0.05:
            value = value -1
        if coeff >= 0.05 and coeff <= 0.15:
            value = value - 0.5
        if coeff >= 0.15 and coeff <= 0.8:
            value = value + 0.5
        if coeff >= 0.8 and coeff <= 1:
            value = value + 1

    if value >= 0 and value < 5:
        coeff = rnd.random()
        if coeff >= 0 and coeff < 0.05:
            value = value -1
        if coeff >= 0.05 and coeff <= 0.25:
            value = value - 0.5
        if coeff >= 0.25 and coeff <= 0.9:
            value = value + 0.5
        if coeff >= 0.9 and coeff <= 1:
            value = value + 1

    if value >= 5 and value < 10:
        coeff = rnd.random()
        if coeff >= 0 and coeff < 0.1:
            value = value -1
        if coeff >= 0.1 and coeff <= 0.4:
            value = value - 0.5
        if coeff >= 0.4 and coeff <= 0.9:
            value = value + 0.5
        if coeff >= 0.9 and coeff <= 1:
            value = value + 1

    if value >= 10 and value < 15:
        coeff = rnd.random()
        if coeff >= 0 and coeff < 0.1:
            value = value -1
        if coeff >= 0.1 and coeff <= 0.55:
            value = value - 0.5
        if coeff >= 0.55 and coeff <= 0.92:
            value = value + 0.5
        if coeff >= 0.92 and coeff <= 1:
            value = value + 1
    
    if value >= 15 and value < 20:
        coeff = rnd.random()
        if coeff >= 0 and coeff < 0.15:
            value = value -1
        if coeff >= 0.15 and coeff <= 0.65:
            value = value - 0.5
        if coeff >= 0.65 and coeff <= 0.95:
            value = value + 0.5
        if coeff >= 0.95 and coeff <= 1:
            value = value + 1

    if value >= 20 and value < 25:
        coeff = rnd.random()
        if coeff >= 0 and coeff < 0.15:
            value = value -1
        if coeff >= 0.15 and coeff <= 0.67:
            value = value - 0.5
        if coeff >= 0.67 and coeff <= 0.95:
            value = value + 0.5
        if coeff >= 0.95 and coeff <= 1:
            value = value + 1

    if value >= 25 and value < 30:
        coeff = rnd.random()
        if coeff >= 0 and coeff < 0.15:
            value = value -1
        if coeff >= 0.15 and coeff <= 0.75:
            value = value - 0.5
        if coeff >= 0.75 and coeff <= 0.95:
            value = value + 0.5
        if coeff >= 0.95 and coeff <= 1:
            value = value + 1

    if value > 30:
        coeff = rnd.random()
        if coeff >= 0 and coeff < 0.17:
            value = value -1
        if coeff >= 0.17 and coeff <= 0.85:
            value = value - 0.5
        if coeff >= 0.85 and coeff <= 0.97:
            value = value + 0.5
        if coeff >= 0.97 and coeff <= 1:
            value = value + 1

    return value


"""     aux_diff_PM10 = Diff_cte1 * (((0 - aUPM10_cell[nCell]) / delta_dist) - ((flux_tang_PM10 * tf_lf) / delta_dist)) * aLong[nEdge]
                aux_diff_PM25 = Diff_cte2 * (((0 - aUPM25_cell[nCell]) / delta_dist) - ((flux_tang_PM25 * tf_lf) / delta_dist)) * aLong[nEdge]

                # FLUJO ADVECTIVO
                aux_adv_PM10 = (G_f + np.absolute(G_f)) * 0.5 * aUPM10_cell[nCell]
                aux_adv_PM25 = (G_f + np.absolute(G_f)) * 0.5 * aUPM25_cell[nCell] """