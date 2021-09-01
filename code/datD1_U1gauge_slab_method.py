import numpy as np

# importo il file analisi_errori.py comune a tutti i moduli
import os
import sys
from inspect import getfile, currentframe
current_dir = os.path.dirname(os.path.abspath(getfile(currentframe())))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from analisi_errori import *

gT = 4  # rapporto fra costante di accoppiamento teoria e temperatura

# carico tutte le sottocartelle cartella
filenames = os.listdir('./DATA/U1gauge_slab')

# definisco un po' di liste utili
dim_tot = []
chi_tot = []
sigma_chi_tot = []
tempo_autocorrelazione_tot = []

for name in filenames:
    # per ogni file trovo la dimensione della simulazione
    dim_tot.append(int(name[3:]))

dim_tot = sorted(dim_tot)

# per ogni dimensione
for dim in dim_tot:
    nx = 20
    nt = dim
    chi = []
    dchi = []
    tempo_autocorrelazione = []

    # come prima cosa controllo che Q**2 = 0 su tutto il reticolo
    simulpath = 'DATA/U1gauge_slab/dim' + str(dim) + '/x0.dat'
    with open(simulpath, 'r') as file:
        Q = np.array([float(sushi) for sushi in next(file).split()])
    print('dim =', dim, '\t sum Q**2 =', np.sum(Q**2))

    # ora studio i sottovolumi
    for x in range(1, 10):
        simulpath = 'DATA/U1gauge_slab/dim' + \
            str(dim) + '/x' + str(x) + '.dat'

        # carico la carica topologica da file
        with open(simulpath, 'r') as file:
            Q = np.array([float(sushi) for sushi in next(file).split()])

        # mi calcolo la suscettibilità
        chi.append((np.mean(Q**2) - np.mean(Q)**2) / gT**2)

        # il tempo di autocorrelazione
        measures = len(Q)
        try:
            tempo_autocorrelazione.append(autocorr_time_definizione(Q))
        except ZeroDivisionError:
            print('ZeroDivisionError')
            tempo_autocorrelazione.append(10**6)

        # e l'errore sulla suscettibilità
        sigma_chi_bootstrap = []
        for k in range(15):
            sigma = sigma_bootstrap_mean(Q, block_size=2 ** k, dimension=0)
            if sigma != 0:
                sigma_chi_bootstrap.append(sigma / gT ** 2)
        dchi.append(sigma_chi_bootstrap)

    chi_tot.append(chi)
    sigma_chi_tot.append(dchi)
    tempo_autocorrelazione_tot.append(tempo_autocorrelazione)

# salvo la simulazione su file
simul_path = 'DATA/U1gauge_slab_results.dat'
with open(simul_path, 'w') as file:
    for dim in dim_tot:
        file.write("%d " % dim)
    file.write("\n")

    for chi in chi_tot:
        for x in chi:
            file.write("%f " % x)
        file.write("\n")
    file.write("\n")

    for tau in tempo_autocorrelazione_tot:
        for x in tau:
            file.write("%f " % x)
        file.write("\n")
    file.write("\n")

    for dchi in sigma_chi_tot:
        for k in dchi:
            for x in k:
                file.write("%f " % x)
            file.write("\n")
        file.write("\n")
