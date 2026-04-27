import numpy as np
import matplotlib.pyplot as plt

def Heaviside(x):
    if x ==0:
        return 0.5
    
    return 1 * (x > 0)

t_min, t_max = -1, 3

time = np.linspace(t_min, t_max, 500)
func_exp = [np.exp(-t) for t in time]
func_exp_max = np.max(func_exp)
func_hs = [np.exp(-t)*Heaviside(t) for t in time]

fig = plt.figure(figsize=(8,10), tight_layout=True)
ax1 = fig.add_subplot(211)
ax1.set_title(r'$E$($t$)')
ax1.set_xlim(t_min, t_max)
ax1.set_ylim(-0.1, 2.0)
ax1.set_xlabel(r'$t$ /s')
ax1.plot(time, func_exp)
ax1.hlines([0], -1, 3, 'black', lw=1, ls='-')
ax1.vlines([0], -0.1, 2, 'black', lw=1, ls='-')
ax1.grid()

ax2 = fig.add_subplot(212)
ax2.set_title(r'$E$($t$) $h_{{s}}$($t$)')
ax2.set_xlim(t_min, t_max)
ax2.set_ylim(-0.1, 2.0)
ax2.set_xlabel(r'$t$ /s')
ax2.plot(time, func_hs)
ax2.hlines([0], -1, 3, 'black', lw=1, ls='-')
ax2.vlines([0], -0.1, 2, 'black', lw=1, ls='-')
ax2.grid()

savefile = './png/Heaviside.png'
fig.savefig(savefile, dpi=300)

plt.show()