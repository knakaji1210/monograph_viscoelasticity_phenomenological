import numpy as np
import matplotlib.pyplot as plt

tauM = 10
tau_min, tau_max = -3, 4

tau = np.logspace(tau_min, tau_max, 100)
log_tau = [np.log10(t) for t in tau]
approxRelax = 2*(tauM/tau)**3 / (1 + (tauM/tau)**2)**2

fig = plt.figure(tight_layout=True)
ax = fig.add_subplot(111)
ax.set_xlim(tau_min, tau_max)
ax.set_ylim(-0.1,1.1)
ax.set_xlabel(r'log($\tau$ /s)')
ax.set_ylabel(r'$H(\tau) / E_i$')
ax.plot(log_tau, approxRelax)
ax.vlines([np.log10(tauM)], 0, 1.1, 'r', ls='--')
ax.grid()
fig.text(0.2,0.8, r'$\tau_M$ = {0:.1f} s'.format(tauM))

savefile = './png/Alfrey_approx_Maxwell_relaxSpectrum_losMod.png'
fig.savefig(savefile, dpi=300)

plt.show()