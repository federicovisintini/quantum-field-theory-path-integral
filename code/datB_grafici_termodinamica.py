import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# importo il file analisi_errori.py comune a tutti i moduli
import os
import sys
from inspect import getfile, currentframe
current_dir = os.path.dirname(os.path.abspath(getfile(currentframe())))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from analisi_errori import *

# Grafici termodinamica campo scalare, energia a massa fissa, RETICOLO ANISOTROPO
mass = 0.05
alpha = 1
size = 14
figsize = (6.3, 5)
save_fig = False

# carico i dati da file
simul_path = 'DATA/camposcalare/massa_fissa_anisotr.dat'
with open(simul_path, 'r') as file:
    nt_arr = np.array([int(nt) for nt in next(file).split()])
    oss = np.array([float(xene) for xene in next(file).split()])
    doss = np.array([float(xene) for xene in next(file).split()])

temp = 1 / nt_arr / mass
min_temp = np.argmin(temp)

ene = (oss - oss[min_temp]) * nt_arr ** 2
dene = doss * nt_arr ** 2
xline = np.linspace(0, max(temp), 1000)

temp, ene, dene = zip(*sorted(zip(temp, ene, dene)))
temp = np.array(temp)
ene = np.array(ene)
dene = np.array(dene)

# mostro l'energia del campo, normalizzata con la temperatura
fig, ax = plt.subplots(1, 1, figsize=figsize)
ax.set_title(r'Energia campo scalare con $\hat{m} = 0.05$', size=size)
ax.errorbar(temp, ene, dene, marker='.', ls='--', alpha=alpha)
ax.plot(xline, np.pi / 6 * np.ones(1000))
ax.set_ylabel(
    r'$\langle \epsilon/T^2 \rangle -  \langle \epsilon/\bar{T}^2 \rangle_{\bar{T}=0}$', size=size)
ax.set_xlabel("T/m", size=size)
ax.tick_params(labelsize=size)
fig.tight_layout()
if save_fig is True:
    fig.savefig("figure/massa_fissa.png")

# Grafici termodinamica campo scalare, energia a massa fissa, ANOMALIA DI TRACCIA


def integrate(x, y, dy):
    val = 0.0
    err = 0.0
    area = [0.0]
    darea = [0.0]
    for i in range(len(x) - 1):
        delta_x = x[i + 1] - x[i]
        val += delta_x * (y[i + 1] + y[i]) / 2.0
        dz = delta_x * (dy[i + 1] + dy[i])
        err += dz ** 2  # sommo errori in quadratura
        area.append(val)
        darea.append(err)

    return(np.array(area), np.sqrt(darea))


# carico i dati da file
simul_path = 'DATA/camposcalare/massa_fissa_trace_anomaly.dat'
with open(simul_path, 'r') as file:
    nt_arr = np.array([int(nt) for nt in next(file).split()])
    O1 = np.array([float(xene) for xene in next(file).split()])
    dO1 = np.array([float(xene) for xene in next(file).split()])


temp = 1 / nt_arr / mass
min_temp = np.argmin(temp)

tr = nt_arr ** 2 * (O1 - O1[min_temp])
dtr = nt_arr ** 2 * dO1
xline = np.linspace(0, max(temp), 1000)
# tr = O1

# ordino l'array delle temperature
temp, tr, dtr = zip(*sorted(zip(temp, tr, dtr)))
temp = np.array(temp)
tr = np.array(tr)
dtr = np.array(dtr)

# calcolo la pressione
pressure, dpressure = integrate(temp, tr / temp, dtr / temp)

# mostro l'anomalia di traccia, normalizzata con la temperatura
fig, ax = plt.subplots(1, 1, figsize=(8 / 1.2, 6 / 1.2))
ax.set_title(r'Energia campo scalare con $\hat{m} = 0.05$', size=size)
ax.errorbar(temp, tr, dtr, marker='.', ls='--',
            label=r'$\epsilon-p$', alpha=alpha)

# mostro la pressione
plt.errorbar(temp, pressure, dpressure, marker='.',
             ls='--', label=r'$p$', alpha=alpha)

# mostro la somma di pressione e anomalia di traccia (energia interna)
plt.errorbar(temp, pressure + tr, dpressure +
             dtr, marker='.', ls='--', label=r'$\epsilon$', alpha=alpha)
plt.errorbar(temp, ene, dene, marker='.', ls='--',
             label=r'$\epsilon$ ret. anisotr.', alpha=alpha)
plt.plot(xline, np.pi / 6 * np.ones(1000))

ax.set_ylabel(
    r'$\langle \epsilon/T^2 \rangle -  \langle \epsilon/\bar{T}^2 \rangle_{\bar{T}=0}$', size=size)
ax.set_xlabel("T/m", size=size)
ax.tick_params(labelsize=size)
ax.legend(frameon=False, fontsize='x-large')

fig.tight_layout()
if save_fig is True:
    fig.savefig("figure/trace_anomaly.png")


# Grafici termodinamica campo scalare, energia a temperatura fissa, RETICOLO ANISOTROPO

def quadr(x, a, b):
    return a + b * x**2


# carico i dati da file
simul_path = 'DATA/camposcalare/massa_variabile_anisotr.dat'
with open(simul_path, 'r') as file:
    nt_arr = np.array([int(nt) for nt in next(file).split()])
    ene2 = np.array([float(xene) for xene in next(file).split()])
    dene2 = np.array([float(xene) for xene in next(file).split()])

    # ordino gli array delle temperature
    nt_arr, ene2, dene2 = zip(*sorted(zip(nt_arr, ene2, dene2)))
    nt_arr = np.array(nt_arr)
    ene2 = np.array(ene2)
    dene2 = np.array(dene2)

masses = 0.1 / nt_arr

ene2 = ene2 * nt_arr ** 2
dene2 = dene2 * nt_arr ** 2
xline = np.linspace(0, max(masses), 1000)

# fit, al variare della massa
opt, cov = curve_fit(quadr, masses[2:], ene2[2:], [0.5, 20], dene2[2:])
print("fit result = {0:.4f} +/- {1:.4f}".format(opt[0], np.sqrt(cov[0, 0])))
print("pi / 6     = {0:.4f}".format(np.pi / 6))

# mostro la costante per T grande al variare della massa
fig, ax = plt.subplots(1, 1, figsize=(8 / 1.2, 6 / 1.2))
ax.set_title(r'Limite continuo a $T/m = 10$', size=size)
ax.errorbar(masses, ene2, dene2, marker='.', ls='')
plt.plot(xline, quadr(xline, *opt))
plt.plot(xline, np.pi / 6 * np.ones(1000))
ax.set_ylabel(r'$\epsilon$', size=size)
ax.set_xlabel(r'$\hat{m}$', size=size)
ax.tick_params(labelsize=size)
fig.tight_layout()
if save_fig is True:
    fig.savefig("figure/massa_variabile.png")

plt.show()
