import numpy as np
import random
from numba import jit


@jit
def force_x(space, time, ix, it, nx, nt):
    force_x = np.exp(
        1j * (time[ix, it] + space[ix, (it + 1) % nt] - time[(ix + 1) % nx, it]))
    force_x += np.exp(1j * (- time[ix, it - 1] +
                            space[ix, it - 1] + time[(ix + 1) % nx, it - 1]))

    return force_x


@jit
def force_t(space, time, ix, it, nx, nt):
    force_t = np.exp(
        1j * (space[ix, it] + time[(ix + 1) % nx, it] - space[ix, (it + 1) % nt]))
    force_t += np.exp(1j * (- space[ix - 1, it] +
                            time[ix - 1, it] + space[ix - 1, (it + 1) % nt]))
    return force_t


@jit
def norm(phi):
    while phi < -np.pi:
        phi = phi + 2 * np.pi
    while phi > np.pi:
        phi = phi - 2 * np.pi
    return phi


@jit
def update_metropolis(space, time, delta, beta, nx, nt):
    acc = 0
    for ix in range(nx):
        for it in range(nt):
            # aggiorno lungo x
            forza_x = force_x(space, time, ix, it, nx, nt)

            phi = space[ix, it]
            phi0 = np.angle(forza_x)
            f = np.absolute(forza_x)
            phi_p = delta * (2 * random.random() - 1) + phi
            r = np.exp(beta * f * (np.cos(phi_p - phi0) - np.cos(phi - phi0)))
            if random.random() < r:
                space[ix, it] = norm(phi_p)
                acc += 1

            # aggiorno lungo t
            forza_t = force_t(space, time, ix, it, nx, nt)

            phi = time[ix, it]
            phi0 = np.angle(forza_t)
            f = np.absolute(forza_t)
            phi_p = delta * (2 * random.random() - 1) + phi
            r = np.exp(beta * f * (np.cos(phi_p - phi0) - np.cos(phi - phi0)))
            if random.random() < r:
                time[ix, it] = norm(phi_p)
                acc += 1

    return space, time, acc / (2 * nx * nt)


@jit
def update_overrelax(space, time, nx, nt):
    for ix in range(nx):
        for it in range(nt):
            # aggiorno lungo x
            forza_x = force_x(space, time, ix, it, nx, nt)

            phi = space[ix, it]
            phi0 = np.angle(forza_x)
            f = np.absolute(forza_x)

            space[ix, it] = norm(2.0 * phi0 - phi)

            # aggiorno lungo t
            forza_t = force_t(space, time, ix, it, nx, nt)

            phi = time[ix, it]
            phi0 = np.angle(forza_t)
            f = np.absolute(forza_t)

            time[ix, it] = norm(2.0 * phi0 - phi)

    return space, time


@jit
def suscett(space, time, nx, nt):
    Q = 0.0
    for ix in range(nx):
        for it in range(nt):
            placquett12 = space[ix, it] + time[(ix + 1) % nx, it]
            placquett12 += -space[ix, (it + 1) % nt] - time[ix, it]
            while placquett12 < - np.pi:
                placquett12 += 2 * np.pi
            while placquett12 > np.pi:
                placquett12 -= 2 * np.pi
            Q += placquett12
    return Q / (2 * np.pi)


@jit
def sotto_suscett(space, time, nx, nt, nx_max, nt_max):
    Q = 0.0
    for ix in range(nx):
        for it in range(nt):
            placquett12 = space[ix, it] + time[(ix + 1) % nx_max, it]
            placquett12 += -space[ix, (it + 1) % nt_max] - time[ix, it]
            while placquett12 < - np.pi:
                placquett12 += 2 * np.pi
            while placquett12 > np.pi:
                placquett12 -= 2 * np.pi
            Q += placquett12
    return Q / (2 * np.pi)
