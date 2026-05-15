# ordinary differential equation of SLS1 model (sinusoidal strain)
# freqSweep version

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

'''
テキストの式(4.11)をベースに組み立てる
'''

def SLS1_sinuStrain(s, t, eamp, af, insMod, k, tau):
# e: strain, s: stress, infMod: relaxation modulus, tau: relaxation time
# ここではeampとafを指定し、この中でeの関数を作り振動歪みを実現
    e = eamp*np.sin(af*t)
    dedt = eamp*af*np.cos(af*t)
    dsdt = (insMod*e + insMod*tau*dedt - k*s)/tau   # (4.11')
    return dsdt

def getNearestIndex2value(list,value):
    index = np.abs(np.array(list) -value).argsort()[0].tolist()
    return index

# variables
try:
    E1 = float(input('modulus 1 [MPa] (default = 1.0 MPa): '))*10**6
except ValueError:
    E1 = 10**6                    # [Pa] modulus
try:
    E2 = float(input('modulus 2 [MPa] (default = 0.5 MPa): '))*10**6
except ValueError:
    E2 = 5*10**5                  # [Pa] modulus
try:
    eta = float(input('viscosity [kPa s] (default = 1000.0 kPa s): '))*10**3
except ValueError:
    eta = 10**6                   # [Pa s] viscosity
tau = eta/E2                  # [s] retardation time
insMod = E1                   # [Pa] instantaneous modulus
infMod = E1*E2/(E1+E2)        # [Pa] equilibrium modulus
k = insMod/infMod

# external sinusoidal strain
try:
    log_af_min = float(input('log(af_min) for forced oscillation (default=-2.0): '))
except ValueError:
    log_af_min = -2.0
try:
    log_af_max = float(input('log(af_max) for forced oscillation (default=1.0): '))
except ValueError:
    log_af_max = 1.0
try:
    num_freq = int(input('number of anglar frequency (default=30): '))
except ValueError:
    num_freq = 30
try:
    eamp = float(input('amplitude for sinusoidal strain [] (default=0.2): '))
except ValueError:
    eamp = 0.2

af_list = np.logspace(log_af_min, log_af_max, num_freq)
aft_list = af_list*tau

dt = 0.01                  # [s] interval time
s0 = 0                      # [] initial strain

e = []              # 周波数掃引の全ての入力信号（歪み）を格納
s = []              # 周波数掃引の全ての出力信号（応力）を格納
samp = []           # 各周波数での出力振幅の最大値を格納
spha = []           # 各周波数での出力信号の位相を格納

for af in af_list:
    tmax = 10*np.pi/af          # [s] duration time
    t = np.arange(0, tmax, dt)    # time during sinusoidal strain
    e_af = eamp*np.sin(af*t)    # 入力信号
    e_af_stat = e_af[int(0.8*len(e_af)):]
    e.extend(e_af)
    # solution of ODE
    sol = odeint(SLS1_sinuStrain, s0, t, args=(eamp,af,insMod,k,tau))
    s_af = sol[:,0]             # [MPa] stress
    s_af_stat = s_af[int(0.8*len(s_af)):]
    s_max = np.max(s_af_stat)
    samp.append(s_max)        
    ind = getNearestIndex2value(s_af_stat,0)    # 出力信号が0になるindexを抽出
    s_pha = (180/np.pi)*np.arcsin(np.abs(e_af_stat[ind])/eamp)
    spha.append(s_pha)  

e = np.array(e)
s = np.array(s)
samp = np.array(samp)
spha = np.array(spha)

# scaling for figure
s = s/10**6                   # MPaスケール
samp = samp/10**6             # MPaスケール

fig = plt.figure(figsize=(8,5), tight_layout=True)
ax1 = fig.add_subplot(111)
ax1.grid()
title_text = "SLS1 model: sinusoidal strain"
fig.suptitle(title_text)
ax1.set_axisbelow(True)
ax1.set_xscale('log')
ax1.set_xlabel(r'$\omega_f$ $\tau$')
ax1.set_ylabel(r'$\sigma_{{amp}}$ /MPa')
ax1.set_ylim(0,1.2*np.max(samp))
ax2 = ax1.twinx()
ax2.grid(ls='dotted')
ax2.set_ylim(0,50)
ax2.set_ylabel(r'$\theta$ /$\degree$')

var_text = r'$\epsilon_{{amp}}$ = {0:.2f}, $E_1$ = {1:.1f} MPa, $E_2$ = {2:.1f} MPa, $\eta$ = {3:.1f} kPa s'.format(eamp,E1/10**6,E2/10**6,eta/10**3)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'd$\sigma$/d$t$ = ($E_i$$\epsilon$ + $E_i$$\tau$ d$\epsilon$/d$t$ - $k$$\sigma$)/$\tau$'
ax1.text(0.1, 0.8, eq_text, transform=ax2.transAxes)
res_text = r'$\tau$ = {0:.1f} s'.format(tau)
ax1.text(0.1, 0.7, res_text, transform=ax2.transAxes)

ax1.plot(aft_list,samp, 'ro-', label=r'$\sigma_{{amp}}$')
ax2.plot(aft_list,spha, 'bo-', label=r'$\theta$')

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2)

savefile = "./png/SLS1_sinuStrain_freqSweep_(epsilon={0:.2f},tau={1:.1f}s,E1={2:.1f}MPa,E2={3:.1f}MPa,eta={4:.1f}kPas).png".format(eamp,tau,E1/10**6,E2/10**6,eta/10**3)
fig.savefig(savefile, dpi=300)

plt.show()