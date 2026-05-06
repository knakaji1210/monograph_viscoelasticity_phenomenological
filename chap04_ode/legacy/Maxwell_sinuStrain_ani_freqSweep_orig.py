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
tau = eta/E
l = 0.1                     # [m] equilibrium length
w = 0.5                     # ratio of dashpot width

# external sinusoidal strain
try:
    log_af_min = float(input('log(af_min) for forced oscillation (default=-1.0): '))
except ValueError:
    log_af_min = -1.0
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

dt = 0.1                   # [s] interval time
s0 = 0                     # [MPa] initial stress

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
    es_af = s_af/E                                   # strain on spring
    ed_af = integral_s/eta                           # strain on dashpot
    es.extend(es_af)
    ed.extend(ed_af)  
    ind = getNearestIndex2value(s_af_stat,0)          # 出力信号が0になるindexを抽出
    s_pha = (180/np.pi)*np.arcsin(np.abs(e_af_stat[ind])/eamp)
    spha.append(s_pha)                             # 出力位相を格納
    spha_ani.extend([n*s_pha for n in ones_list])   # 最大振幅を格納（アニメーション用） 

t_ani = np.arange(0, i_ani*dt, dt)
e = np.array(e)
s = np.array(s)
es = np.array(es)
ed = np.array(ed)

# 直列なのでl0=2lとなっている
el = e*2*l                                  # [m] elongation
el_s = es*2*l                              # [m] elongation of spring
el_d = ed*2*l                              # [m] elongation of dashpot

# scaling for figure
s = s/10**6                     # MPaスケール

fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111, xlabel='$t$ /s')
ax.grid()
title_text = "Maxwell model: sinusoidal strain"
ax.set_title(title_text)
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
ax.set_xlim(-0.05,0.4)
ax.set_ylim(-5,5)
# for common
y_0 = [0, 0]
ax.plot([0, 0.08*l],y_0, c='b')
ax.plot(0,0,'ro', markersize='10')
ax.plot([0,2*l],[-3,-3], c='g')
ax.plot([0,0],[-2.8,-3.2], c='g')
ax.plot([2*l,2*l],[-2.8,-3.2], c='g')
# for dashpot
x_d1 = [0.08*l, 0.92*l]
y_d1 = [w, w]
y_d2 = [-w, -w]
ax.plot(x_d1,y_d1, c='b')
ax.plot(x_d1,y_d2, c='b')
ax.plot([0.08*l,0.08*l],[w,-w], c='b')
rect = patches.Rectangle(xy=(0.08*l, -w), width=0.83*l, height=2*w, facecolor='y')
ax.add_patch(rect)

var_text = r'$\epsilon_{{amp}}$ = {0:.2f}, $E$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s'.format(eamp,E/10**6,eta/10**3)
ax.text(0.4, 0.9, var_text, transform=ax.transAxes)
eq_text = r'd$\sigma$/d$t$ = -$\sigma$/$\tau$ + $E$d$\epsilon$/d$t$'
ax.text(0.4, 0.8, eq_text, transform=ax.transAxes)
res_text = r'$\tau$ = {0:.1f} s'.format(tau)
ax.text(0.4, 0.7, res_text, transform=ax.transAxes)
ax.text(0.35, 0.15, '$l_0$', transform=ax.transAxes)
ax.text(0.75, 0.28, '$\epsilon$ (input)', transform=ax.transAxes)

# for common
rod, = ax.plot([],[], 'b', animated=True)
point, = ax.plot([],[], 'ro', markersize='10', animated=True)
# for spring
rod_sp, = ax.plot([],[], 'b', animated=True)
triangle, = ax.plot([],[], 'b', animated=True)
# for dashpot
rod_da, = ax.plot([],[], 'b', animated=True)
damper, = ax.plot([],[], 'b', lw=4, animated=True)
# ここでは[],[]としているが、下で***.set_dataで実際の値を入れている

# for strain
strain, = ax.plot([],[], 'b', lw=2, animated=True)
arrow_p, = ax.plot([],[], 'b', marker=9, markersize='10', animated=True)
arrow_n, = ax.plot([],[], 'b', marker=8, markersize='10', animated=True)

samp_template = r'$\sigma_{{amp}}$ = %.2f MPa'
samp_text = ax.text(0.75, 0.52, '', transform=ax.transAxes)
spha_template = r'$\theta$ = %.1f$\degree$'
spha_text = ax.text(0.75, 0.45, '', transform=ax.transAxes)
af_template = r'$\omega_f$ = %.2f s$^{{-1}}$'
af_text = ax.text(0.1, 0.8, '', transform=ax.transAxes)

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return rod, point, triangle, rod_da, damper, strain, arrow_p, arrow_n, samp_text, spha_text, af_text, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    # for common
    x_o = 7*l/4 + el[i]
    x_rod = [x_o, l/4+x_o]
    rod.set_data(x_rod,y_0)
    point.set_data([l/4+x_o],[0])
    # for spring
    x_sp = x_o-el_s[i]-l/2
    x_tri = np.linspace(x_sp, x_o,100)
    y_tri = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri - x_sp)/(x_o-x_sp)-np.pi/2+0.1))-1)
    triangle.set_data(x_tri,y_tri)
    # for dashpot
    x_da = l/2 + el_d[i]
    x_rod_da = [x_da, x_da+3*l/4]
    x_damp = x_da
    y_damp = 0.7*w
    x_damper = [x_damp, x_damp]
    y_damper = [y_damp, -y_damp]
    rod_da.set_data(x_rod_da,y_0)
    damper.set_data(x_damper,y_damper)
    # for strain
    x_strain = [2*l, 2*l + el[i]]
    strain.set_data(x_strain,[-2,-2])
    if e[i] > 0:
        arrow_p.set_data([2*l + el[i]],[-2])
    else:
        arrow_n.set_data([2*l + el[i]],[-2])
    samp_text.set_text(samp_template % samp_ani[i])
    spha_text.set_text(spha_template % spha_ani[i])
    af_text.set_text(af_template % af_ani[i])
    # for time
    time_text.set_text(time_template % (i*dt))
    return rod, point, triangle, rod_da, damper, strain, arrow_p, arrow_n, samp_text, spha_text, af_text, time_text

'''
y_triの中の重要部分は
x_tri1 = np.linspace(a, b,100)
のとき
(xtri - a)/(b - a)
になる 
'''

f = np.arange(0, i_ani)
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f, 
                    init_func=init, blit=True, interval=frame_int, repeat=False)

savefile = "./gif/Maxwell_sinuStrain_freqSweep_(epsilon={0:.2f},tau={1:.1f}s,mod={2:.1f}MPa,eta={3:.1f}kPas).gif".format(eamp,tau,E/10**6,eta/10**3)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()