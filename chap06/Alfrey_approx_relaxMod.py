# Alfrey approximation for relaxation modulus

import numpy as np
import matplotlib.pyplot as plt

time = 10
tau_min, tau_max = 10**(-3), 40

tau = np.linspace(tau_min, tau_max, 200)
func = np.exp(-time/tau)

fig = plt.figure(tight_layout=True)
ax = fig.add_subplot(111)
ax.set_xlim(0, tau_max)
ax.set_ylim(-0.1,1.1)
ax.set_xlabel(r'$\tau$ /s')
ax.set_ylabel(r'$\exp(-t/\tau)$')
ax.plot(tau, func)
ax.hlines([0], tau_min, time, 'r', ls='--')
ax.hlines([1], time, tau_max, 'r', ls='--')
ax.hlines([1/np.exp(1)], tau_min, time, 'g', ls='--')
ax.vlines([time], 0, 1, 'r', ls='--')
ax.grid()
fig.text(0.22,0.86, r'$t$ = 10 s')
fig.text(0.24,0.22, r'$\tau < t$')
fig.text(0.35,0.22, r'$t < \tau$')
fig.text(0.15,0.42, r'1/e')

savefile = './png/Alfrey_approx_relaxMod.png'
fig.savefig(savefile,dpi=300)

plt.show()