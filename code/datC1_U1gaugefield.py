import numpy as np

# importo il file analisi_errori.py comune a tutti i moduli
import os
import sys
from inspect import getfile, currentframe
current_dir = os.path.dirname(os.path.abspath(getfile(currentframe())))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from analisi_errori import *

# rapporto fra costante di accoppiamento teoria e temperatura
gT = 4

# carico tutti i files nella cartella
filenames = os.listdir('./DATA/U1gaugefield')

# definisco un po' di liste utili
dim_tot = []
chi_tot = []
sigma_chi_tot = []
tempo_autocorrelazione_tot = []

for name in filenames:
    # per ogni file trovo la dimensione della simulazione
    dim = int(name[3:-4])
    dim_tot.append(dim)

    # mi calcolo il beta usato (no 1/T), a g^2 costante
    beta = dim ** 2 / gT ** 2

    # carico la carica topologica da file
    with open('DATA/U1gaugefield/' + name, 'r') as file:
        Q = np.array([int(sushi) for sushi in next(file).split()])

    # calcolo il tempo di autocorr di Q, la suscett, e l'errore (no correl)
    measures = len(Q)
    try:
        tempo_autocorrelazione = autocorr_time_definizione(Q)
    except ZeroDivisionError:
        print('ZeroDivisionError')
        tempo_autocorrelazione = 10 ** 6

    tempo_autocorrelazione_tot.append(tempo_autocorrelazione)

    # calcolo la suscetticità
    chi_tot.append((np.mean(Q**2) - np.mean(Q)**2) / gT**2)

    # calcolo l'errore sulla sucettibilità con bootstrap migliorato
    sigma_chi_bootstrap = []

    for k in range(15):
        sigma = sigma_bootstrap_mean(Q, block_size=2 ** k, dimension=0)
        if sigma != 0:
            sigma_chi_bootstrap.append(sigma / gT ** 2)
    sigma_chi_tot.append(sigma_chi_bootstrap)


# salvo la simulazione su file
simul_path = 'DATA/U1gaugefield_results.dat'
with open(simul_path, 'w') as file:
    for dim in dim_tot:
        file.write("%d " % dim)
    file.write("\n")

    for chi in chi_tot:
        file.write("%f " % chi)
    file.write("\n")

    for tau in tempo_autocorrelazione_tot:
        file.write("%f " % tau)
    file.write("\n")

    for dchi in sigma_chi_tot:
        for k in dchi:
            file.write("%f " % k)
        file.write("\n")
