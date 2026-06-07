# Alfrey approximation for complex modulus

import numpy as np
import matplotlib.pyplot as plt

angFreq = 0.1
tau_min, tau_max = 0, 40

tau = np.linspace(tau_min, tau_max, 200)

try:
    select = int(input('Selection (storage modulus: 0, loss modulus: 1): '))
except ValueError:
    select = 0

if select == 0:
    func = angFreq**2 * tau**2 / (1 + angFreq**2 * tau**2)
    ylabel = r'$\omega^2\tau^2 / (1 + \omega^2\tau^2)$'
    savefile = './png/Alfrey_approx_strMod.png'
else:
    func = 1 / (1 + angFreq**2 * tau**2)
    ylabel = r'$1 / (1 + \omega^2\tau^2)$'
    savefile = './png/Alfrey_approx_losMod.png'

fig = plt.figure(tight_layout=True)
ax = fig.add_subplot(111)
ax.set_xlim(tau_min, tau_max)
ax.set_ylim(-0.1,1.1)
ax.set_xlabel(r'$\tau$ /s')
ax.set_ylabel(ylabel)
ax.plot(tau, func)
ax.hlines([0], tau_min, 1/angFreq, 'r', ls='--')
ax.hlines([1], 1/angFreq, tau_max, 'r', ls='--')
ax.hlines([1/2], tau_min, 1/angFreq, 'g', ls='--')
ax.vlines([1/angFreq], 0, 1, 'r', ls='--')
ax.grid()
fig.text(0.18,0.86, r'1/$\omega$ = 10 s')
fig.text(0.22,0.22, r'$\tau < 1/\omega$')
fig.text(0.35,0.22, r'$1/\omega < \tau$')
fig.text(0.15,0.56, '0.5')

fig.savefig(savefile,dpi=300)

plt.show()