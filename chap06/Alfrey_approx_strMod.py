import numpy as np
import matplotlib.pyplot as plt

freq = 0.1
tau_min, tau_max = 0, 40

tau = np.linspace(tau_min, tau_max, 200)
func = [freq**2 * t**2/(1 + freq**2 * t**2) for t in tau]

fig = plt.figure(tight_layout=True)
ax = fig.add_subplot(111)
ax.set_xlim(tau_min, tau_max)
ax.set_ylim(-0.1,1.1)
ax.set_xlabel('tau /s')
ax.set_ylabel('freq^2 * tau^2/(1 + freq^2 * tau^2)')
ax.plot(tau, func)
ax.hlines([0], tau_min, 1/freq, 'r', ls='--')
ax.hlines([1], 1/freq, tau_max, 'r', ls='--')
ax.hlines([1/2], tau_min, 1/freq, 'g', ls='--')
ax.vlines([1/freq], 0, 1, 'r', ls='--')
ax.grid()
fig.text(0.35,0.2, '1/freq = 10 s')
fig.text(0.075,0.48, '0.5')

fig.savefig('./png/approx_form2.png')

plt.show()