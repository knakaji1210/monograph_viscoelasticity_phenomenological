# ordinary differential equation of spring (sinusoidal strain)
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

# external sinusoidal strain
try:
    eamp = float(input('amplitude for sinusoidal strain [] (default=0.25): '))
except ValueError:
    eamp = 0.25
try:
    T = float(input('peirod for sinusoidal strain [s] (default=2.0): '))
except ValueError:
    T = 2.0

af = 2*np.pi/T

tmax = 5*T                       # [s] duration time
dt = 0.01                       # [s] interval time
t_a = np.arange(0, tmax, dt)    # time after step stress
t_b = np.arange(-0.1*tmax,0,dt) # time before step stress
t = np.concatenate([t_b,t_a])   # whole time 
zeros = np.zeros(len(t_b))
e_a = np.array([eamp*np.sin(af*t) for t in t_a])
e = np.concatenate([zeros,e_a]) # whole strain

# solution of ODE（ここではそれをする必要はない）
s = E*e

# scaling for figure
s = s/10**6                     # MPaスケール

fig = plt.figure(figsize=(8,10), tight_layout=True)
title_text = "spring (Hooke's elasticity): sinusoidal strain"
fig.suptitle(title_text)
ax1 = fig.add_subplot(211, xlabel='$t$ /s')
ax2 = fig.add_subplot(212, xlabel='$t$ /s')
ax1.grid()
ax2.grid()
ax1.set_axisbelow(True)
ax2.set_axisbelow(True)
ax1.set_ylabel('strain, $\epsilon$')
ax2.set_ylabel('stress, $\sigma$ /MPa')
ax1.set_ylim(-1.5*np.max(e),1.5*np.max(e))
ax2.set_ylim(-1.5*np.max(s),1.5*np.max(s))

var_text = r'$\epsilon_{{amp}}$ = {0:.2f}, $T$ = {1:.1f} s, $E$ = {2:.1f} MPa'.format(eamp,T,E/10**6)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'$\sigma$ = $E\epsilon$'
ax2.text(0.1, 0.9, eq_text, transform=ax2.transAxes)

ax1.plot(t, e, 'b', label='$\epsilon$ (input)')
ax1.legend(loc='upper right')
ax2.plot(t, s, 'r', label='$\sigma$ (output)')
ax2.legend(loc='upper right')

savefile = "./png/spring_sinuStrain_(epsilon={0:.2f},T={1:.1f}s,mod={2:.1f}MPa).png".format(eamp,T,E/10**6)
fig.savefig(savefile, dpi=300)

plt.show()