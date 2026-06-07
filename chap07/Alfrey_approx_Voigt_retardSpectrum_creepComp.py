# Alfrey approximation example of retardation spectrum of Voigt model by creep compliance

import numpy as np
import matplotlib.pyplot as plt

tauV = 10
tau_min, tau_max = -3, 4

tau = np.logspace(tau_min, tau_max, 100)
log_tau = np.log10(tau)
approxRelax = tau * np.exp(-tau/tauV) / tauV

fig = plt.figure(tight_layout=True)
ax = fig.add_subplot(111)
ax.set_xlim(tau_min, tau_max)
ax.set_ylim(-0.1,1.1)
ax.set_xlabel(r'log($\tau$ /s)')
ax.set_ylabel(r'$L(\tau) E_\infty$')
ax.plot(log_tau, approxRelax)
ax.vlines([np.log10(tauV)], 0, 1.1, 'r', ls='--')
ax.grid()
fig.text(0.2,0.8, r'$\tau_V$ = {0:.1f} s'.format(tauV))

savefile = './png/Alfrey_approx_Voigt_retardSpectrum_creepComp.png'
fig.savefig(savefile, dpi=300)

plt.show()