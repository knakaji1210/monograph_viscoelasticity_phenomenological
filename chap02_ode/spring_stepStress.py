# ordinary differential equation of spring (step stress)
'''
バネ要素単独では常微分方程式は不要
'''

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# variables
try:
    E = float(input('modulus [MPa] (default = 0.2 MPa): '))*10**6
except ValueError:
    E = 2*10**5             # [Pa] modulus

# initial condition
try:
    s0 = float(input('step stress [MPa] (default = 0.05 MPa): '))*10**6
except ValueError:
    s0 = 0.05*10**6         # [Pa] step stress
e0 = 0                      # [] initial strain

tmax = 10                   # [s] duration time
dt = 0.01                   # [s] interval time
t_a = np.arange(0, tmax, dt)    # time after step stress
t_b = np.arange(-0.1*tmax,0,dt) # time before step stress
t = np.concatenate([t_b,t_a])   # whole time 
zeros = np.zeros(len(t_b))
ones = np.ones(len(t_a))
s = np.concatenate([zeros,ones*s0])

# solution of ODE（ここではそれをする必要はない）
e = s/E

# scaling for figure
s = s/10**6                     # MPaスケール

fig = plt.figure(figsize=(8,10), tight_layout=True)
title_text = "spring (Hooke's elasticity): step stress"
fig.suptitle(title_text)
ax1 = fig.add_subplot(211, xlabel='$t$ /s')
ax2 = fig.add_subplot(212, xlabel='$t$ /s')
ax1.grid()
ax2.grid()
ax1.set_axisbelow(True)
ax2.set_axisbelow(True)
ax1.set_ylabel('stress, $\sigma$ /MPa')
ax2.set_ylabel('strain, $\epsilon$')
ax1.set_ylim(-0.1*np.max(s),1.5*np.max(s))
ax2.set_ylim(-0.1*np.max(e),1.5*np.max(e))

var_text = r'$\sigma_0$ = {0:.2f} MPa, $E$ = {1:.1f} MPa'.format(s0/10**6,E/10**6)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'$\epsilon$ = $\sigma_0/E$'
ax2.text(0.1, 0.9, eq_text, transform=ax2.transAxes)

ax1.plot(t, s, 'r', label='$\sigma$ (input)')
ax1.legend(loc='upper right')
ax2.plot(t, e, 'b', label='$\epsilon$ (output)')
ax2.legend(loc='upper right')

savefile = "./png/spring_stepStress.png"
fig.savefig(savefile, dpi=300)

plt.show()