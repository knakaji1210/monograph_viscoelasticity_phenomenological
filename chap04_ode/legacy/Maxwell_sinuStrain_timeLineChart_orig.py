# ordinary differential equation of Maxwell model (sinusoidal strain) with animation
# timeChart version

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.animation import FuncAnimation

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

def getNearestIndex2value(list,value):
    index = np.abs(np.array(list) -value).argsort()[0].tolist()
    return index

# variables
try:
    E = float(input('modulus [MPa] (default = 0.2 MPa): '))*10**6
except ValueError:
    E = 2*10**5                 # [Pa] modulus
try:
    eta = float(input('viscosity [kPa s] (default = 100.0 kPa s): '))*10**3
except ValueError:
    eta = 10**5                 # [Pa s] viscosity
tau = eta/E

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
    tmax = 10*np.pi/af         # [s] duration time
    t = np.arange(0, tmax, dt) # time after step stress
    i_ani += len(t)
    zeros = np.zeros(len(t))
    ones = np.ones(len(t))
    ones_list = ones.tolist()  # アニメーション用（各周波数でのアニメーション期間中に同じ数値をずっと表示させるため）
    af_ani.extend([n*af for n in ones_list])  # 入力周波数を格納（アニメーション用）
    e_af = eamp*np.sin(af*t)                  # 入力信号
    e_af_stat = e_af[int(0.8*len(e_af)):]             # 後半部分を抽出（前半は過渡応答を含むから）
    e.extend(e_af)                                    # 入力信号（アニメーション用）
    sol = odeint(Maxwell_sinuStrain, s0, t, args=(eamp,af,E,tau))  # ODEの解を求めている
    s_af = sol[:, 0]                                  # [s]が出てくる
    s_af_stat = s_af[int(0.8*len(s_af)):]             # 後半部分を抽出（前半は過渡応答を含むから）
    s_max = np.max(s_af_stat)                         # 後半部分の最大値（最大振幅と見做す）
    samp.append(s_max/10**6)                          # 最大振幅を格納
    samp_ani.extend([n*s_max/10**6 for n in ones_list])     # 最大振幅を格納（アニメーション用）
    s.extend(s_af)                                    # 出力信号（アニメーション用）
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

fig = plt.figure(figsize=(8,5), tight_layout=True)
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
ax1.set_ylim(-3*eamp, 3*eamp)
ax2.set_ylim(-2*samp_max, 2*samp_max)
ax1.grid(False)
ax2.grid(True)
ax1.set_axisbelow(True)
ax2.set_axisbelow(True)
ax1.set_xlabel('$t$ /s')
ax1.set_ylabel('strain, $\epsilon$')
ax2.set_ylabel('stress, $\sigma$ /MPa')

var_text = r'$\epsilon_{{amp}}$ = {0:.2f}, $E$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s'.format(eamp,E/10**6,eta/10**3)
ax1.text(0.1, 0.25, var_text, transform=ax1.transAxes)
eq_text = r'd$\sigma$/d$t$ = -$\sigma$/$\tau$ + $E$d$\epsilon$/d$t$'
ax1.text(0.1, 0.15, eq_text, transform=ax1.transAxes)
res_text = r'$\tau$ = {0:.1f} s'.format(tau)
ax1.text(0.1, 0.05, res_text, transform=ax1.transAxes)

input, = ax1.plot([], [], 'b', animated=True, label='$\epsilon$ (input)')
output, = ax2.plot([], [], 'r', animated=True, label='$\sigma$ (output)')
# ここでは[],[]としているが、下で***.set_data([0, l + x[i]], [0, 0])で実際の値を入れている

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2)

samp_template = r'$\sigma_{{amp}}$ = %.2f MPa'
samp_text = ax1.text(0.35, 0.9, '', transform=ax1.transAxes)
spha_template = r'$\theta$ = %.1f$\degree$'
spha_text = ax1.text(0.35, 0.8, '', transform=ax1.transAxes)
af_template = r'$\omega_f$ = %.2f s$^{{-1}}$'
af_text = ax1.text(0.1, 0.8, '', transform=ax1.transAxes)

time_template = '$t$ = %.1f s'
time_text = ax1.text(0.1, 0.9, '', transform=ax1.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return input, output, samp_text, spha_text, af_text, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    if i < 200:
        ax1.set_xlim(0,t_ani[200])
        input.set_data([t_ani[:i]],[e[:i]])
        output.set_data([t_ani[:i]],[s[:i]])
    else:
        ax1.set_xlim(t_ani[i-200],t_ani[i])
        input.set_data([t_ani[i-200:i]],[e[i-200:i]])
        output.set_data([t_ani[i-200:i]],[s[i-200:i]])
    samp_text.set_text(samp_template % samp_ani[i])
    spha_text.set_text(spha_template % spha_ani[i])
    af_text.set_text(af_template % af_ani[i])
    # for time
    time_text.set_text(time_template % (i*dt))
    return input, output, samp_text, spha_text, af_text, time_text

f = np.arange(0, i_ani)
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f, 
                    init_func=init, blit=True, interval=frame_int, repeat=False)

savefile = "./gif/Maxwell_sinuStrain_timeChart_(tau={0:.1f}s).gif".format(tau)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()