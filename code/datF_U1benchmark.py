import matplotlib.pyplot as plt
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
dim = 8
clr = 'C0'
beta = dim ** 2 / gT ** 2

# carico la carica topologica da file
with open('DATA/U1gaugefield_benchmark.dat', 'r') as file:
    Q = np.array([int(sushi) for sushi in next(file).split()])

# calcolo il tempo di autocorr di Q, la suscett, e l'errore (no correl)
measures = len(Q)
tempo_autocorrelazione = autocorr_time_definizione(Q)

chi = (np.mean(Q**2) - np.mean(Q)**2) / gT**2
chi_previsto = 1 / (4 * np.pi ** 2)  # risultato teorico a T = 0

# calcolo l'errore sulla sucettibilità con bootstrap migliorato
sigma_chi_bootstrap = []
for k in range(13):
    sigma_chi_bootstrap.append(sigma_bootstrap_mean(
        Q, block_size=2 ** k, dimension=0) / gT**2)

# e prendo il massimo dell'errore
sigma_chi_bootstrap = np.array(sigma_chi_bootstrap)
sigma_chi = np.mean(sigma_chi_bootstrap[8:])

plt.figure('carica_topologica_dim' + str(dim))
plt.plot(range(measures), Q, marker='.', linestyle='', color=clr)

plt.figure('bootstrap_migliorato' + str(dim))
plt.plot(range(13), sigma_chi_bootstrap,
         marker='.', linestyle='', color='C0')
plt.hlines(sigma_chi, 0, 12, colors=['C1'])
plt.plot(range(8, 13), sigma_chi_bootstrap[8:],
         marker='.', linestyle='', color='C1')
plt.title('bootstrap su blocchi di grandezza $2^k$')
plt.xlabel('k: blocks of size $2^k$')
plt.ylabel(r'$\sigma_k$')

print('chi = {0:.6f} +/- {1:.6f}'.format(chi, sigma_chi))
print('chi previsto =', chi_previsto)
print('tempo autocorr {0:.4f}'.format(tempo_autocorrelazione))

plt.show()
