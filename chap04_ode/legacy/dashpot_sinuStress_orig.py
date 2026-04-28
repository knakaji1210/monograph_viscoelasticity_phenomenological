# ordinary differential equation of dashpot (sinusoidal stress)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

'''
テキストの式(2.5)をベースに組み立てる
'''

def dashpot_sinuStress(e, t, samp, af, eta):
# e: strain, s: stress, eta: viscosity
# ここではsampとafを指定し、この中でsの関数を作り振動応力を実現
    s = samp*np.sin(af*t)
    dedt = s/eta        # (1.5)
    return dedt

# variables
try:
    eta = float(input('viscosity [kPa s] (default = 500.0 kPa s): '))*10**3
except ValueError:
    eta = 5*10**5             # [Pa s] viscosity

# external sinusoidal strain
try:
    samp = float(input('amplitude for sinusoidal stress [MPa] (default=0.2): '))*10**6
except ValueError:
    samp = 0.2*10**6
try:
    T = float(input('peirod for sinusoidal strain [s] (default=2.0): '))
except ValueError:
    T = 2.0

af = 2*np.pi/T
e0 = 0.0

tmax = 5*T                  # [s] duration time
dt = T/250                   # [s] interval time
t_a = np.arange(0, tmax, dt)    # time after step stress
t_b = np.arange(-0.1*tmax,0,dt) # time before step stress
t = np.concatenate([t_b,t_a])   # whole time 
zeros = np.zeros(len(t_b))
s_a = np.array([samp*np.sin(af*t) for t in t_a])
s = np.concatenate([zeros,s_a]) # whole stress

# solution of ODE
sol = odeint(dashpot_sinuStress, e0, t_a, args=(samp,af,eta))
e = np.concatenate([zeros,sol[:,0]])

# scaling for figure
s = s/10**6                     # MPaスケール

fig = plt.figure(figsize=(8,10), tight_layout=True)
title_text = "dashpot (Newton's viscosity): sinusoidal stress"
fig.suptitle(title_text)
ax1 = fig.add_subplot(211, xlabel='$t$ /s')
ax2 = fig.add_subplot(212, xlabel='$t$ /s')
ax1.grid()
ax2.grid()
ax1.set_axisbelow(True)
ax2.set_axisbelow(True)
ax1.set_ylabel('stress, $\sigma$ /MPa')
ax2.set_ylabel('strain, $\epsilon$')
ax1.set_ylim(-1.5*np.max(s),1.5*np.max(s))
ax2.set_ylim(-1.5*np.max(e),1.5*np.max(e))

var_text = r'$\sigma_{{amp}}$ = {0:.2f} MPa, $T$ = {1:.1f} s, $\eta$ = {2:.1f} kPa s'.format(samp/10**6,T,eta/10**3)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'd$\epsilon$/d$t$ = $\sigma$/$\eta$'
ax2.text(0.1, 0.9, eq_text, transform=ax2.transAxes)

ax1.plot(t, s, 'r', label='$\sigma$ (input)')
ax1.legend(loc='upper right')
ax2.plot(t, e, 'b', label='$\epsilon$ (output)')
ax2.legend(loc='upper right')

savefile = "./png/dashpot_sinuStress_(sigma={0:.2f}MPa,T={1:.1f}s,eta={2:.1f}kPas).png".format(samp,T,eta/10**3)
fig.savefig(savefile, dpi=300)

plt.show()