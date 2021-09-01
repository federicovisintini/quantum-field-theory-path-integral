import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def line(x, a, b):
    return a * x + b


save_fig = False
size = 14  # figure

# carico i dati da file
simul_path = 'DATA/U1gaugefield_results.dat'
with open(simul_path, 'r') as file:
    dim = np.array([int(x) for x in next(file).split()])
    chi = np.array([float(x) for x in next(file).split()])
    tau = np.array([float(x) for x in next(file).split()])
    num = len(chi)
    err = np.ndarray(num, dtype=list)
    for i in range(num):
        err[i] = np.array([float(delta) for delta in next(file).split()])

# riordino gli array
dim, chi, tau, err = zip(*sorted(zip(dim, chi, tau, err)))
dim = np.array(dim)
chi = np.array(chi)
tau = np.array(tau)
err = np.array(err)

# risultato teorico a T = 0
chi_previsto = 1 / (4 * np.pi ** 2)

# calcolo l'errore con bootstrap migliorato
sigma_chi = np.ndarray(num)
for i in range(num):
    sigma_chi[i] = np.mean(err[i][5 + i:15])

# faccio la figura della suscettività e dell'autocorrelazione
fig, [ax_susc, ax_corr] = plt.subplots(
    2, 1, figsize=(6.3, 9.5), sharex=True)

ax_susc.set_title(r'$\chi$ for standard metropolis', size=size)
ax_susc.set_ylabel(r"$\chi$", size=size)
ax_susc.tick_params(labelsize=size)
ax_susc.errorbar(dim**2, chi, sigma_chi, marker='.', linestyle='', color='C0')
ax_susc.hlines(chi_previsto, min(dim)**2, max(dim)
               ** 2, color='C1', linestyle='--')

ax_corr.set_title('Critical slowing down', size=size)
ax_corr.set_ylabel(r'$\tau_{int}$', size=size)
ax_corr.set_xlabel(r'$N_s \times N_t$', size=size)
ax_corr.tick_params(labelsize=size)
ax_corr.set_yscale('log')
ax_corr.scatter(dim**2, tau, color='C0')

popt, pcov = curve_fit(line, dim**2, np.log10(tau))
ax_corr.plot(dim**2, 10 ** (dim ** 2 *
                            popt[0] + popt[1]), color='gray', linestyle='--')

fig.tight_layout()
if save_fig is True:
    fig.savefig("figure/suscettivity.png")

# Bootstrap migliorato
dim_boot = 8
index = np.where(dim == dim_boot)[0][0]

fig, ax_bsm = plt.subplots(figsize=(6.3, 5))
ax_bsm.scatter(range(15), err[index][:], color='C0')
ax_bsm.hlines(y=sigma_chi[index], xmin=0., xmax=14, color='C1', linestyle='--')
ax_bsm.scatter(range(dim_boot, 15), err[index][dim_boot:15], color='C1')

ax_bsm.set_title("Bootstrap migliorato", size=size)
ax_bsm.set_ylabel(r"$\sigma_k$", size=size)
ax_bsm.set_xlabel("$k$", size=size)
ax_bsm.tick_params(labelsize=size)
ax_bsm.set_yscale('log')

fig.tight_layout()
if save_fig is True:
    fig.savefig("figure/bootstrap_migliorato8.png")


plt.show()
