import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def line(x, a, b):
    return a * x + b


def rainbow(x, chi_vero):
    return x * (1 - x) * chi_vero


save_fig = False
size = 14  # figure

# carico i dati da file
simul_path = 'DATA/U1gauge_slab_results.dat'
with open(simul_path, 'r') as file:
    dim = np.array([int(x) for x in next(file).split()])

    chi = np.ndarray(len(dim), dtype=list)
    for i in range(len(dim)):
        chi[i] = np.array([float(x) for x in next(file).split()])
    next(file)

    tau = np.ndarray([len(dim), 9])
    for i in range(len(dim)):
        tmp = np.array([float(x) for x in next(file).split()])
        for x in range(len(tmp)):
            tau[i, x] = tmp[x]
    next(file)

    err = np.ndarray([len(dim), 9], dtype=list)
    for i in range(len(dim)):
        for subvol in range(9):
            err[i, subvol] = np.array([float(x) for x in next(file).split()])
        next(file)

nx = 20
nt = dim

# risultato teorico a T = 0
chi_previsto = 1 / (4 * np.pi ** 2)

# calcolo l'errore con bootstrap migliorato
sigma_chi = np.ndarray([len(dim), 9])
for i in range(len(dim)):
    for x in range(9):
        sigma_chi[i, x] = np.mean(err[i, x][8: 15])

# Bootstrap migliorato X = 0.1
dim_boot = 17
x = 0

index = np.where(dim == dim_boot)[0][0]

fig, ax_bsm = plt.subplots(figsize=(6.3, 5))
ax_bsm.scatter(range(15), err[index, x][:], color='C0')
ax_bsm.hlines(y=sigma_chi[index, x], xmin=0.,
              xmax=15, color='C1', linestyle='--')
ax_bsm.scatter(range(8, 15), err[index, x][8:], color='C1')

ax_bsm.set_title("Bootstrap migliorato, $x = 0.1$", size=size)
ax_bsm.set_ylabel(r"$\sigma_k$", size=size)
ax_bsm.set_xlabel("$k$", size=size)
ax_bsm.tick_params(labelsize=size)
ax_bsm.set_yscale('log')

fig.tight_layout()
if save_fig is True:
    fig.savefig("figure/bootstrap_migliorato_slab17.png")

# fit suscettibilità
x = np.linspace(0.1, 0.9, 9)
x_line = np.linspace(0, 1, 1000)
chi_fit = []
dchi_fit = []

for i in range(len(dim)):
    popt, pcov = curve_fit(rainbow, x, chi[i], [0.02533], sigma_chi[i])
    chi_fit.append(popt)
    dchi_fit.append(np.sqrt(pcov[0, 0]))

print('coeff dim: {0}, chi = {1:.5f} +/- {2:.5f}'.format(
    dim[index], chi_fit[index][0], dchi_fit[index]))
fig, ax_fit = plt.subplots(figsize=(6.3, 5))

ax_fit.errorbar(x, chi[index], sigma_chi[index],
                marker='.', linestyle='', color='black')
ax_fit.plot(x_line, rainbow(x_line, chi_fit[index]), color='red')
ax_fit.set_title(r'Fit di $\chi$', size=size)
ax_fit.set_xlabel("x", size=size)
ax_fit.set_ylabel(r"$\chi_s$", size=size)
ax_fit.tick_params(labelsize=size)

fig.tight_layout()
if save_fig is True:
    fig.savefig("figure/suscettivity_fit17.png")

fig, [ax_susc, ax_corr] = plt.subplots(
    2, 1, figsize=(6.3, 9.5), sharex=True)

ax_susc.set_title(r'$\chi$ for slab method, $T \sim 0$')
ax_susc.set_ylabel(r"$\chi$", size=size)
ax_susc.tick_params(labelsize=size)
ax_susc.errorbar(nx * nt, chi_fit, dchi_fit,
                 marker='.', linestyle='', color='C0')
ax_susc.hlines(chi_previsto, min(nt) * nx, max(nt)
               * nx, color='C1', linestyle='--')

ax_corr.set_title('Autocorrelation time, slab method')
ax_corr.set_ylabel(r'$\tau_{int}$', size=size)
ax_corr.set_xlabel(r'$N_s \times N_t$', size=size)
ax_corr.tick_params(labelsize=size)
ax_corr.set_yscale('log')
ax_corr.scatter(nx * nt, tau[:, 0], label='x=0.1')
ax_corr.scatter(nx * nt, tau[:, 1], label='x=0.2')
ax_corr.legend(frameon=False, fontsize='x-large')

popt, pcov = curve_fit(line, np.log10(nx * nt), np.log10(tau[:, 0]))
ax_corr.plot(nx * nt, (nx * nt) **
             popt[0] * 10 ** popt[1], color='C0', linestyle='--')

print('coeff x=0.1: {0:.2f} +/- {1:.2f}'.format(popt[0], np.sqrt(pcov[0, 0])))
popt, pcov = curve_fit(line, np.log10(nx * nt), np.log10(tau[:, 1]))
ax_corr.plot(nx * nt, (nx * nt) **
             popt[0] * 10 ** popt[1], color='C1', linestyle='--')

print('coeff x=0.2: {0:.2f} +/- {1:.2f}'.format(popt[0], np.sqrt(pcov[0, 0])))
fig.tight_layout()

if save_fig is True:
    fig.savefig("figure/suscettivity_slab.png")

plt.show()
