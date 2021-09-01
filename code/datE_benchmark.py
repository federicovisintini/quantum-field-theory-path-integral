from scalar2d import *
import numpy as np

# importo il file analisi_errori.py comune a tutti i moduli
import os
import sys
from inspect import getfile, currentframe
current_dir = os.path.dirname(os.path.abspath(getfile(currentframe())))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from analisi_errori import *

# Analisi termodinamica campo scalare, energia a massa fissa

mass = 0.05
nt = 30
simul_path = 'DATA/camposcalare/benchmark.dat'
print('Analisi a massa fissa - benchmark')

# lista dove metto le osservabili a meno di dividere per T**2
oss = []
doss = []
tr = []
dtr = []


xene_mass = []
xene_spat = []
xene_temp = []

# carico dati da file
with open(simul_path, 'r') as file:
    xene_mass += [float(xene) for xene in next(file).split()]
    xene_spat += [float(xene) for xene in next(file).split()]
    xene_temp += [float(xene) for xene in next(file).split()]

xene_mass = np.array(xene_mass)
xene_spat = np.array(xene_spat)
xene_temp = np.array(xene_temp)

measures = len(xene_mass)

tr.append(xene_mass.mean())
dtr.append(xene_mass.std() / np.sqrt(measures) *
           np.sqrt(1 + 2 * autocorr_time_definizione(xene_mass)))

# Valore delle osservabili
# print("<m^2 phi^2>             = {0:.6f} +/- {1:.6f}".format(O1, dO1))
# print("<(phi(n+x) - phi(n))^2> = {0:.6f} +/- {1:.6f}".format(O2, dO2))
# print("<(phi(n+t) - phi(n))^2> = {0:.6f} +/- {1:.6f}".format(O3, dO3))

# calcolo l'energia (poi la normalizzerò per T^2)
energy = (xene_mass + xene_spat - xene_temp) / 2
denergy = energy.std() / np.sqrt(measures) * np.sqrt(1 +
                                                     2 * autocorr_time_definizione(energy))

oss.append(energy.mean())
doss.append(denergy)

# print("Nt = {2}, energ = {0:.6f} +/- {1:.6}".format(energy.mean(), denergy, nt))
# print("Nt = {2}, trace = {0:.6f} +/- {1:.6}".format(tr[-1], dtr[-1], nt))

print(denergy * 10 ** 5)
print(dtr[-1] * 10 ** 5)
