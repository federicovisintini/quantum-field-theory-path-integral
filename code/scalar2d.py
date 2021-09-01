import numpy as np
import random
from numba import jit


@jit
def energy(field, mass2):
    nx, nt = field.shape
    nvol = nx * nt
    xene_mass = 0.0
    xene_spat = 0.0
    xene_temp = 0.0
    for ix in range(nx):
        for it in range(nt):
            phi = field[ix, it]
            force_s = field[(ix + 1) % nx, it]
            force_t = field[ix, (it + 1) % nt]

            xene_mass += mass2 * phi**2
            xene_spat += - 2.0 * phi * force_s + 2.0 * phi**2
            xene_temp += - 2.0 * phi * force_t + 2.0 * phi**2

    return xene_mass / float(nvol), xene_spat / float(nvol), xene_temp / float(nvol)


@jit
def update_heatbath(field, mass2):
    mass2p4 = mass2 + 4.0
    nx, nt = field.shape
    for ix in range(nx):
        for it in range(nt):
            force = 0.0
            phi = field[ix, it]
            force += field[(ix + 1) % nx, it]
            force += field[ix - 1, it]
            force += field[ix, (it + 1) % nt]
            force += field[ix, it - 1]

            # variance of the gaussian distribution is 1/(m^2 + 4)
            sigma2 = 1.0 / mass2p4

            # average of the gaussian distribution is force/(m^2 + 4)
            aver = force * sigma2

            # BOX MULLER ALGORITHM
            x = np.sqrt(sigma2) * np.sqrt(-2.0 * np.log(random.random()))
            y = x * np.cos(2.0 * np.pi * random.random()) + aver

            # update field
            field[ix, it] = y

    return field


@jit
def update_overrelax(field, mass2):
    mass2p4 = mass2 + 4.0
    nx, nt = field.shape
    for ix in range(nx):
        for it in range(nt):
            force = 0.0
            phi = field[ix, it]
            force += field[(ix + 1) % nx, it]
            force += field[ix - 1, it]
            force += field[ix, (it + 1) % nt]
            force += field[ix, it - 1]

            # average of the gaussian distribution is force / (m^2 + 4)
            aver = force / mass2p4

            field[ix, it] = 2.0 * aver - phi

    return field
