# ordinary differential equation of Maxwell model (sinusoidal strain) with animation
# freqSweep version

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.animation import FuncAnimation

'''
テキストの式(1.11)をベースに組み立てる
'''

def Maxwell_sinuStrain(s, t, eamp, af, E, tau):
# e: strain, s: stress, E: modulus, tau: relaxation time
# ここではeampとafを指定し、この中でeの関数を作り振動歪みを実現
    e = eamp*np.sin(af*t)
    dedt = eamp*af*np.cos(af*t)
    dsdt = E*dedt-s/tau         # (1.11)
    return dsdt

def getNearestIndex2value(list,value):
    index = np.abs(np.array(list) -value).argsort()[0].tolist()
    return index

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

dt = 0.01                   # [s] interval time
s0 = 0                     # [MPa] initial stress

e = []              # 周波数掃引の全ての入力信号（歪み）を格納
s = []              # 周波数掃引の全ての出力信号（応力）を格納
es = []             # 周波数掃引の全ての出力信号（バネの歪み）を格納
ed = []             # 周波数掃引の全ての出力信号（ダッシュポットの歪み）を格納
samp = []           # 各周波数での出力振幅の最大値を格納
spha = []           # 各周波数での出力信号の位相を格納

for af in af_list:
    tmax = 10*np.pi/af         # [s] duration time
    t = np.arange(0, tmax, dt) # time during sinusoidal strain
    e_af = eamp*np.sin(af*t)   # 入力信号
    e_af_stat = e_af[int(0.8*len(e_af)):]             # 後半部分を抽出（前半は過渡応答を含むから）
    e.extend(e_af)                                    # 入力信号（アニメーション用）
    sol = odeint(Maxwell_sinuStrain, s0, t, args=(eamp,af,E,tau))  # ODEの解を求めている
    s_af = sol[:, 0]                                  # [s]が出てくる
    s_af_stat = s_af[int(0.8*len(s_af)):]             # 後半部分を抽出（前半は過渡応答を含むから）
    s_max = np.max(s_af_stat)                         # 後半部分の最大値（最大振幅と見做す）
    samp.append(s_max/10**6)                          # 最大振幅を格納
    integral_s = np.array([s_af[:k+1].sum()*dt for k in range(len(s_af))])     # 簡易的なsの積分
    es_af = s_af/E                                   # strain on spring
    ed_af = integral_s/eta                           # strain on dashpot
    es.extend(es_af)
    ed.extend(ed_af)  
    ind = getNearestIndex2value(s_af_stat,0)          # 出力信号が0になるindexを抽出
    s_pha = (180/np.pi)*np.arcsin(np.abs(e_af_stat[ind])/eamp)
    spha.append(s_pha)                             # 出力位相を格納

e = np.array(e)
s = np.array(s)
#es = np.array(es)
#ed = np.array(ed)

# scaling for figure
s = s/10**6                     # MPaスケール

samp_max = np.max(samp)

fig = plt.figure(figsize=(8,5), tight_layout=True)
ax1 = fig.add_subplot(111)
ax1.grid()
title_text = "Maxwell model: sinusoidal strain"
ax1.set_title(title_text)
ax1.set_axisbelow(True)
ax1.set_xlabel(r'$\omega_f$ $\tau$')
ax1.set_ylim(0,1.2*samp_max)
ax1.set_ylabel(r'$\sigma_{{amp}}$ /MPa')
ax2 = ax1.twinx()
ax2.grid(ls='dotted')
ax2.set_ylim(-10,100)
ax2.set_ylabel(r'$\theta$ /$\degree$')

var_text = r'$\epsilon_{{amp}}$ = {0:.2f}, $E$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s'.format(eamp,E/10**6,eta/10**3)
ax1.text(0.05, 0.5, var_text, transform=ax1.transAxes)
eq_text = r'd$\sigma$/d$t$ = -$\sigma$/$\tau$ + $E$d$\epsilon$/d$t$'
ax1.text(0.05, 0.4, eq_text, transform=ax1.transAxes)
res_text = r'$\tau$ = {0:.1f} s'.format(tau)
ax1.text(0.05, 0.3, res_text, transform=ax1.transAxes)

ax1.plot(aft_list,samp, 'ro-', label=r'$\sigma_{{amp}}$')
ax1.set_xscale('log')
ax2.plot(aft_list,spha, 'bo-', label=r'$\theta$')

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2)

savefile = "./png/Maxwell_sinuStrain_freqSweep_(epsilon={0:.2f},tau={1:.1f}s,mod={2:.1f}MPa,eta={3:.1f}kPas).png".format(eamp,tau,E/10**6,eta/10**3)
fig.savefig(savefile, dpi=300)

plt.show()