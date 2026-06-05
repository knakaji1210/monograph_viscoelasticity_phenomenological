import numpy as np
import matplotlib.pyplot as plt

time = 10
tau_min, tau_max = 0, 40

tau = np.linspace(tau_min, tau_max, 200)
func = [np.exp(-time/t) for t in tau]

fig = plt.figure(tight_layout=True)
ax = fig.add_subplot(111)
ax.set_xlim(tau_min, tau_max)
ax.set_ylim(-0.1,1.1)
ax.set_xlabel('tau /s')
ax.set_ylabel('exp(-t/tau)')
ax.plot(tau, func)
ax.hlines([0], tau_min, time, 'r', ls='--')
ax.hlines([1], time, tau_max, 'r', ls='--')
ax.hlines([1/np.exp(1)], tau_min, time, 'g', ls='--')
ax.vlines([time], 0, 1, 'r', ls='--')
ax.grid()
fig.text(0.35,0.2, 't = 10 s')
fig.text(0.075,0.38, '1/e')

fig.savefig('./png/approx_form1.png')

plt.show()