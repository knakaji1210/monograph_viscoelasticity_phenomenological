# ordinary differential equation of dashpot (step stress)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

'''
テキストの式(1.5)をベースに組み立てる
'''

def dashpot_stepStress(e, t, s, eta):
# e: strain, s: stress, eta: viscosity
# ここでは下でargsとしてs=s0を入れてステップ応力を実現
    dedt = s/eta        # (1.5)
    return dedt

# variables
try:
    eta = float(input('viscosity [kPa s] (default = 500.0 kPa s): '))*10**3
except ValueError:
    eta = 5*10**5             # [Pa s] viscosity

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

# solution of ODE
sol = odeint(dashpot_stepStress, e0, t_a, args=(s0,eta))
e = np.concatenate([zeros,sol[:,0]])

# scaling for figure
s = s/10**6                     # MPaスケール

fig = plt.figure(figsize=(8,10), tight_layout=True)
title_text = "dashpot (Newton's viscosity): step stress"
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

var_text = r'$\sigma_0$ = {0:.2f} MPa, $\eta$ = {1:.1f} kPa s'.format(s0/10**6,eta/10**3)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'd$\epsilon$/d$t$ = $\sigma_0$/$\eta$'
ax2.text(0.1, 0.9, eq_text, transform=ax2.transAxes)

ax1.plot(t, s, 'r', label='$\sigma$ (input)')
ax1.legend(loc='upper right')
ax2.plot(t, e, 'b', label='$\epsilon$ (output)')
ax2.legend(loc='upper right')

savefile = "./png/dashpot_stepStress.png"
fig.savefig(savefile, dpi=300)

plt.show()