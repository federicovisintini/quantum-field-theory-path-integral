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
filenames = os.listdir('./DATA/camposcalare/m' + str(mass))
print('Analisi a massa fissa')

# lista dove metto le osservabili a meno di dividere per T**2
oss = []
doss = []
tr = []
dtr = []

temperature = []

for name in filenames:

    nt = name[2:-4]
    temperature.append(int(nt))
    xene_mass = []
    xene_spat = []
    xene_temp = []

    # carico dati da file
    simul_path = 'DATA/camposcalare/m' + str(mass) + '/' + name
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

    print("Nt = {2}, energ = {0:.5f} +/- {1:.5}".format(energy.mean(), denergy, nt))
    print("Nt = {2}, trace = {0:.5f} +/- {1:.5}".format(tr[-1], dtr[-1], nt))

# salvo la simulazione su file
simul_path = 'DATA/camposcalare/massa_fissa_anisotr.dat'
with open(simul_path, 'w') as file:
    for nt in temperature:
        file.write("%d " % nt)
    file.write("\n")

    for xene in oss:
        file.write(str(xene) + " ")
    file.write("\n")

    for dxene in doss:
        file.write(str(dxene) + " ")

# salvo la trace anomaly su file
simul_path = 'DATA/camposcalare/massa_fissa_trace_anomaly.dat'
with open(simul_path, 'w') as file:
    for nt in temperature:
        file.write("%d " % nt)
    file.write("\n")

    for xene in tr:
        file.write(str(xene) + " ")
    file.write("\n")

    for dxene in dtr:
        file.write(str(dxene) + " ")


# Analisi termodinamica campo scalare, energia a temperatura fissa
print('\nAnalisi a temperatura fissa')
filenames = os.listdir('./DATA/camposcalare/T10')

# lista dove metto le osservabili a meno di dividere per T**2
oss = []
doss = []
temperature = []

for name in filenames:

    nt = name[2:-4]
    temperature.append(int(nt))
    # print('Nt =', nt)

    # carico dati da file
    simul_path = 'DATA/camposcalare/T10/' + name
    with open(simul_path, 'r') as file:
        xene_mass = np.array([float(xene) for xene in next(file).split()])
        xene_spat = np.array([float(xene) for xene in next(file).split()])
        xene_temp = np.array([float(xene) for xene in next(file).split()])

    measures = len(xene_mass)

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

    print(
        "Nt = {2}, energy = {0:.5f} +/- {1:.5}".format(energy.mean(), denergy, nt))

# salvo la simulazione su file
simul_path = 'DATA/camposcalare/massa_variabile_anisotr.dat'
with open(simul_path, 'w') as file:
    for nt in temperature:
        file.write("%d " % nt)
    file.write("\n")

    for xene in oss:
        file.write(str(xene) + ' ')
    file.write("\n")

    for dxene in doss:
        file.write(str(dxene) + ' ')
