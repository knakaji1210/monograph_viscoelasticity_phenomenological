# ordinary differential equation of Voigt model (sinusoidal stress)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

'''
テキストの式(1.21)をベースに組み立てる
'''

def Voigt_sinuStress(e, t, samp, af, tau):
# e: strain, s: stress, tau: retardation time
# ここではsampとafを指定し、この中でsの関数を作り振動応力を実現
    s = samp*np.sin(af*t)
    dedt = (s/E - e)/tau        # (1.21)
    return dedt

# variables
try:
    E = float(input('modulus [MPa] (default = 0.2 MPa): '))*10**6
except ValueError:
    E = 2*10**5                 # [Pa] modulus
try:
    eta = float(input('viscosity [kPa s] (default = 500.0 kPa s): '))*10**3
except ValueError:
    eta = 5*10**5               # [Pa s] viscosity
tau = eta/E                     # [s] retardation time

# external sinusoidal strain
try:
    samp = float(input('amplitude for sinusoidal stress [MPa] (default=0.1): '))*10**6
except ValueError:
    samp = 0.1*10**6
try:
    T = float(input('peirod for sinusoidal strain [s] (default=2.0): '))
except ValueError:
    T = 2.0

af = 2*np.pi/T
e0 = 0                      # [] initial strain

tmax = 5*T                  # [s] duration time
dt = 0.01                   # [s] interval time
t_a = np.arange(0, tmax, dt)    # time after step stress
t_b = np.arange(-0.1*tmax,0,dt) # time before step stress
t = np.concatenate([t_b,t_a])   # whole time 
zeros = np.zeros(len(t_b))
s_a = np.array([samp*np.sin(af*t) for t in t_a])
s = np.concatenate([zeros,s_a]) # whole stress

# solution of ODE
sol = odeint(Voigt_sinuStress, e0, t_a, args=(samp,af,tau))
e = np.concatenate([zeros,sol[:,0]])        # [] strain
dedt = np.array([0.0]+[(e[k+1]-e[k])/(t[k+1]-t[k]) for k in range(len(e)-1)])     # 簡易的なde/dt
s_s = E*e                                   # stress on spring
s_d = eta*dedt                              # stress on dashpot

# scaling for figure
s = s/10**6                     # MPaスケール
s_s = s_s/10**6                 # MPaスケール
s_d = s_d/10**6                 # MPaスケール

fig = plt.figure(figsize=(8,10), tight_layout=True)
title_text = "Voigt model: sinusoidal stress"
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

var_text = r'$\sigma_{{amp}}$ = {0:.2f} MPa, $T$ = {1:.1f} s, $E$ = {2:.1f} MPa, $\eta$ = {3:.1f} kPa s'.format(samp/10**6,T,E/10**6,eta/10**3)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'd$\epsilon$/d$t$ = ($\sigma$/$E$ - $\epsilon$)/$\tau$'
ax2.text(0.1, 0.9, eq_text, transform=ax2.transAxes)
res_text = r'$\tau$ = {0:.1f} s'.format(tau)
ax2.text(0.1, 0.8, res_text, transform=ax2.transAxes)

ax1.plot(t, s, 'r', label='$\sigma$ (input)')
ax1.plot(t, s_s, 'g', ls="dashed", label='$\sigma$ (spring)')
ax1.plot(t, s_d, 'y', ls="dashed", label='$\sigma$ (dashpot)')
ax1.legend(loc='upper right')
ax2.plot(t, e, 'b', label='$\epsilon$ (output)')
ax2.legend(loc='upper right')


savefile = "./png/Voigt_sinuStress_(sigma={0:.2f}MPa,T={1:.1f}s,mod={2:.1f}MPa,eta={3:.1f}kPas).png".format(samp/10**6,T,E/10**6,eta/10**3)
fig.savefig(savefile, dpi=300)

plt.show()