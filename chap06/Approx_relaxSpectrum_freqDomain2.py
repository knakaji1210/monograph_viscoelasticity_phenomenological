import numpy as np
import matplotlib.pyplot as plt

T = 10
tau_min, tau_max = -3, 4

tau = np.logspace(tau_min, tau_max, 100)
print(tau)
log_tau = [np.log10(t) for t in tau]
approxRelax = [2*(T/t)**3 / (1 + (T/t)**2)**2 for t in tau]

fig = plt.figure(tight_layout=True)
ax = fig.add_subplot(111)
ax.set_xlim(tau_min, tau_max)
ax.set_ylim(-0.1,1.1)
ax.set_xlabel('log(tau /s)')
ax.set_ylabel('H(tau)')
ax.plot(log_tau, approxRelax)
ax.vlines([np.log10(T)], 0, 1.1, 'r', ls='--')
ax.grid()
fig.text(0.2,0.3, 'tau_M = {0} s'.format(T))

fig.savefig('./png/approx_relax_spectrum_strMod.png')

plt.show()