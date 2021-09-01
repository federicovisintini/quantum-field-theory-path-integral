import numpy as np
import random
from numba import jit


@jit
def energy(field, mass2):
    nx, ny, nz, nt = field.shape
    nvol = nx * ny * nz * nt
    xene_mass = 0.0
    xene_spat = 0.0
    xene_temp = 0.0
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for it in range(nt):
                    phi = field[ix, iy, iz, it]
                    force_s = field[(ix + 1) % nx, iy, iz, it]
                    force_s += field[ix, (iy + 1) % ny, iz, it]
                    force_s += field[ix, iy, (iz + 1) % nz, it]
                    force_t = field[ix, iy, iz, (it + 1) % nt]

                    xene_mass += mass2 * phi**2
                    xene_spat += - 2.0 * phi * force_s + 6.0 * phi**2
                    xene_temp += - 2.0 * phi * force_t + 2.0 * phi**2

    return xene_mass / float(nvol), xene_spat / float(nvol), xene_temp / float(nvol)


@jit
def update_heatbath(field, mass2):
    mass2p8 = mass2 + 8.0
    nx, ny, nz, nt = field.shape
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for it in range(nt):
                    force = 0.0
                    phi = field[ix, iy, iz, it]
                    force += field[(ix + 1) % nx, iy, iz, it]
                    force += field[ix - 1, iy, iz, it]
                    force += field[ix, (iy + 1) % ny, iz, it]
                    force += field[ix, iy - 1, iz, it]
                    force += field[ix, iy, (iz + 1) % nz, it]
                    force += field[ix, iy, iz - 1, it]
                    force += field[ix, iy, iz, (it + 1) % nt]
                    force += field[ix, iy, iz, it - 1]

                    # variance of the gaussian distribution is 1/(m^2 + 8)
                    sigma2 = 1.0 / mass2p8

                    # average of the gaussian distribution is force/(m^2 + 8)
                    aver = force * sigma2

                    # BOX MULLER ALGORITHM
                    x = np.sqrt(sigma2) * np.sqrt(-2.0 *
                                                  np.log(random.random()))
                    y = x * np.cos(2.0 * np.pi * random.random()) + aver

                    # update field
                    field[ix, iy, iz, it] = y

    return field


@jit
def update_overrelax(field, mass2):
    mass2p8 = mass2 + 8
    nx, ny, nz, nt = field.shape
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for it in range(nt):
                    force = 0.0
                    phi = field[ix, iy, iz, it]
                    force += field[(ix + 1) % nx, iy, iz, it]
                    force += field[ix - 1, iy, iz, it]
                    force += field[ix, (iy + 1) % ny, iz, it]
                    force += field[ix, iy - 1, iz, it]
                    force += field[ix, iy, (iz + 1) % nz, it]
                    force += field[ix, iy, iz - 1, it]
                    force += field[ix, iy, iz, (it + 1) % nt]
                    force += field[ix, iy, iz, it - 1]

                    # average of the gaussian distribution is force / (m^2 + 8)
                    aver = force / mass2p8

                    field[ix, iy, iz, it] = 2.0 * aver - phi

    return
