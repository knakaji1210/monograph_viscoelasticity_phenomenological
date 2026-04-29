# ordinary differential equation of Maxwell model (sinusoidal strain)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation

'''
テキストの式(2.11)をベースに組み立てる
'''

def Maxwell_sinuStrain(s, t, eamp, af, E, tau):
# e: strain, s: stress, E: modulus, tau: relaxation time
# ここではeampとafを指定し、この中でeの関数を作り振動歪みを実現
    e = eamp*np.sin(af*t)
    dedt = eamp*af*np.cos(af*t)
    dsdt = E*dedt-s/tau         # (1.11)
    return dsdt

# variables
try:
    E = float(input('modulus [MPa] (default = 0.2 MPa): '))*10**6
except ValueError:
    E = 2*10**5                 # [Pa] modulus
try:
    eta = float(input('viscosity [kPa s] (default = 500.0 kPa s): '))*10**3
except ValueError:
    eta = 5*10**5               # [Pa s] viscosity
tau = eta/E                     # [s] relaxation time

# external sinusoidal strain
try:
    eamp = float(input('amplitude for sinusoidal strain [] (default=0.2): '))
except ValueError:
    eamp = 0.2
try:
    T = float(input('peirod for sinusoidal strain [s] (default=2.0): '))
except ValueError:
    T = 2.0

af = 2*np.pi/T

s0 = 0                      # [] initial strain

tmax = 4*T                  # [s] duration time
dt = 0.01                   # [s] interval time
t_a = np.arange(0, tmax, dt)    # time after step stress
t_b = np.arange(-0.1*tmax,0,dt) # time before step stress
t = np.concatenate([t_b,t_a])   # whole time 
zeros = np.zeros(len(t_b))
e_a = np.array([eamp*np.sin(af*t) for t in t_a])
e = np.concatenate([zeros,e_a]) # whole strain

# solution of ODE
sol = odeint(Maxwell_sinuStrain, s0, t_a, args=(eamp,af,E,tau))
s = np.concatenate([zeros,sol[:,0]])        # [] stress
integral_s = np.array([s[:k+1].sum()*dt for k in range(len(s))])     # 簡易的なsの積分
e_s = s/E                                   # strain on spring
e_d = integral_s/eta                        # strain on dashpot

# scaling for figure
s = s/10**6                     # MPaスケール

fig = plt.figure(figsize=(8,7), tight_layout=True)
title_text = "Maxwell model: sinusoidal strain"
fig.suptitle(title_text)
ax1 = fig.add_subplot(211, xlabel='$t$ /s')
ax2 = fig.add_subplot(212, xlabel='$t$ /s')
ax1.grid()
ax2.grid()
ax1.set_axisbelow(True)
ax2.set_axisbelow(True)
ax1.set_ylabel('strain, $\epsilon$')
ax2.set_ylabel('stress, $\sigma$ /MPa')
ax1.set_xlim(t[0], t[-1])
ax2.set_xlim(t[0], t[-1])
ax1.set_ylim(-1.5*np.max(e),1.5*np.max(e))
ax2.set_ylim(-1.5*np.max(s),1.5*np.max(s))

var_text = r'$\epsilon_{{amp}}$ = {0:.2f}, $T$ = {1:.1f} s, $E$ = {2:.1f} MPa, $\eta$ = {3:.1f} kPa s'.format(eamp,T,E/10**6,eta/10**3)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'd$\sigma$/d$t$ = -$\sigma$/$\tau$ + $E$d$\epsilon$/d$t$'
ax2.text(0.1, 0.9, eq_text, transform=ax2.transAxes)
res_text = r'$\tau$ = {0:.1f} s'.format(tau)
ax2.text(0.1, 0.8, res_text, transform=ax2.transAxes)

line_strain, = ax1.plot([], [], 'b', label='$\epsilon$ (input)')
line_strain_s, = ax1.plot([], [], 'g', ls="dashed", label='$\epsilon$ (spring)')
line_strain_d, = ax1.plot([], [], 'y', ls="dashed", label='$\epsilon$ (dashpot)')
ax1.legend(loc='upper right')
line_stress, =  ax2.plot([], [], 'r', label='$\sigma$ (output)')
ax2.legend(loc='upper right')

# 3. アニメーション更新関数
def animate(i):
    # 応力データの更新
    line_strain.set_data(t[:i], e[:i])
    line_strain_s.set_data(t[:i], e_s[:i])
    line_strain_d.set_data(t[:i], e_d[:i])

    # 歪みデータの更新
    line_stress.set_data(t[:i], s[:i])
    return line_strain, line_strain_s, line_strain_d, line_stress

# アニメーション実行
ani = animation.FuncAnimation(fig, animate, frames=len(t), interval=20, blit=True)

savefile = "./gif/Maxwell_sinuStrain_(T={0:.1f}s).gif".format(T)
ani.save(savefile, writer='pillow', fps=50)

plt.show()