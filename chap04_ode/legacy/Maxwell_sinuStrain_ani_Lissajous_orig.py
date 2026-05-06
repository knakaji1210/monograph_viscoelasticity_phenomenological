# ordinary differential equation of Maxwell model (sinusoidal strain) with animation
# Lissajous version

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
    E = float(input('modulus [MPa] (default = 1.0 MPa): '))*10**6
except ValueError:
    E = 10**6               # [Pa] modulus
try:
    eta = float(input('viscosity [kPa s] (default = 500.0 kPa s): '))*10**3
except ValueError:
    eta = 5*10**5           # [Pa s] viscosity
tau = eta/E                 # [s] relaxation time

# external sinusoidal strain
try:
    log_af_min = float(input('log(af_min) for forced oscillation (default=0.0): '))
except ValueError:
    log_af_min = 0.0
try:
    log_af_max = float(input('log(af_max) for forced oscillation (default=1.0): '))
except ValueError:
    log_af_max = 1.0
try:
    num_freq = int(input('number of anglar frequency (default=5): '))
except ValueError:
    num_freq = 5
try:
    eamp = float(input('amplitude for sinusoidal strain [] (default=0.2): '))
except ValueError:
    eamp = 0.2

af_list = np.logspace(log_af_min, log_af_max, num_freq)

dt = 0.05           # [s] interval time
s0 = 0              # [MPa] initial stress

i_ani = 0
e = []              # 周波数掃引の全ての入力信号（歪み）を格納
s = []              # 周波数掃引の全ての出力信号（応力）を格納
es = []             # 周波数掃引の全ての出力信号（バネの歪み）を格納
ed = []             # 周波数掃引の全ての出力信号（ダッシュポットの歪み）を格納
samp = []           # 各周波数での出力振幅の最大値を格納
spha = []           # 各周波数での出力信号の位相を格納
af_ani = []         # 入力角周波数を格納（アニメーション用）
samp_ani = []       # 出力振幅の最大値を格納（アニメーション用）
spha_ani = []       # 出力信号の位相を格納（アニメーション用）

for af in af_list:
    tmax = 10*np.pi/af                                # [s] duration time
    t = np.arange(0, tmax, dt)                        # time after step stress
    i_ani += len(t)
    zeros = np.zeros(len(t))
    ones = np.ones(len(t))
    ones_list = ones.tolist()       # アニメーション用（各周波数でのアニメーション期間中に同じ数値をずっと表示させるため）
    af_ani.extend([n*af for n in ones_list])          # 入力周波数を格納（アニメーション用）
    e_af = eamp*np.sin(af*t)                          # 入力信号
    e_af_stat = e_af[int(0.8*len(e_af)):]             # 後半部分を抽出（前半は過渡応答を含むから）
    e.extend(e_af)                                    # 入力信号（アニメーション用）
    sol = odeint(Maxwell_sinuStrain, s0, t, args=(eamp,af,E,tau))  # ODEの解を求めている
    s_af = sol[:, 0]                                  # [s]が出てくる
    s_af_stat = s_af[int(0.8*len(s_af)):]             # 後半部分を抽出（前半は過渡応答を含むから）
    s_max = np.max(s_af_stat)                         # 後半部分の最大値（最大振幅と見做す）
    samp.append(s_max/10**6)                          # 最大振幅を格納
    samp_ani.extend([n*s_max/10**6 for n in ones_list])     # 最大振幅を格納（アニメーション用）
    s.extend(s_af)                                          # 出力信号（アニメーション用）
    integral_s = np.array([s_af[:k+1].sum()*dt for k in range(len(s_af))])     # 簡易的なsの積分
    es_af = s_af/E                                    # strain on spring
    ed_af = integral_s/eta                            # strain on dashpot
    es.extend(es_af)
    ed.extend(ed_af)  
    ind = getNearestIndex2value(s_af_stat,0)          # 出力信号が0になるindexを抽出
    s_pha = (180/np.pi)*np.arcsin(np.abs(e_af_stat[ind])/eamp)
    spha.append(s_pha)                                # 出力位相を格納
    spha_ani.extend([n*s_pha for n in ones_list])     # 最大振幅を格納（アニメーション用） 

t_ani = np.arange(0, i_ani*dt, dt)
e = np.array(e)
s = np.array(s)
#es = np.array(es)      # ここでは使っていない
#ed = np.array(ed)      # ここでは使っていない

# scaling for figure
s = s/10**6                     # MPaスケール

samp_max = np.max(samp)

fig = plt.figure(figsize=(8,6), tight_layout=True)
ax = fig.add_subplot(111)
ax.set_xlim(-1.5*eamp, 1.5*eamp)
ax.set_ylim(-1.5*samp_max, 1.5*samp_max)
ax.grid(True)
ax.set_axisbelow(True)
ax.set_xlabel('strain, $\epsilon$')
ax.set_ylabel('stress, $\sigma$ /MPa')

var_text = r'$\epsilon_{{amp}}$ = {0:.2f}, $E$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s'.format(eamp,E/10**6,eta/10**3)
ax.text(0.5, 0.25, var_text, transform=ax.transAxes)
eq_text = r'd$\sigma$/d$t$ = -$\sigma$/$\tau$ + $E$d$\epsilon$/d$t$'
ax.text(0.5, 0.15, eq_text, transform=ax.transAxes)
res_text = r'$\tau$ = {0:.1f} s'.format(tau)
ax.text(0.5, 0.05, res_text, transform=ax.transAxes)

lissajous, = plt.plot([], [], 'r', animated=True)
# ここでは[],[]としているが、下で***.set_data([0, l + x[i]], [0, 0])で実際の値を入れている

samp_template = r'$\sigma_{{amp}}$ = %.2f MPa'
samp_text = ax.text(0.35, 0.9, '', transform=ax.transAxes)
spha_template = r'$\theta$ = %.1f$\degree$'
spha_text = ax.text(0.35, 0.8, '', transform=ax.transAxes)
af_template = r'$\omega_f$ = %.2f s$^{{-1}}$'
af_text = ax.text(0.1, 0.8, '', transform=ax.transAxes)

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return lissajous, samp_text, spha_text, af_text, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    if i < 50:
        lissajous.set_data([e[:i]],[s[:i]])
    else:
        lissajous.set_data([e[i-50:i]],[s[i-50:i]])
    samp_text.set_text(samp_template % samp_ani[i])
    spha_text.set_text(spha_template % spha_ani[i])
    af_text.set_text(af_template % af_ani[i])
    # for time
    time_text.set_text(time_template % (i*dt))
    return lissajous, samp_text, spha_text, af_text, time_text

f = np.arange(0, i_ani)
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f, 
                    init_func=init, blit=True, interval=frame_int, repeat=False)

savefile = "./gif/Maxwell_sinuStrain_Lissajous_(epsilon={0:.2f},tau={1:.1f}s,mod={2:.1f}MPa,eta={3:.1f}kPas).gif".format(eamp,tau,E/10**6,eta/10**3)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()