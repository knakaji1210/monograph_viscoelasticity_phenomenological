# ordinary differential equation of Maxwell model (sinusoidal strain) with animation

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
tau = eta/E                     # [s] relaxation time
l = 0.1                         # [m] equilibrium length
w = 0.5                         # ratio of dashpot width

# external sinusoidal strain
try:
    eamp = float(input('amplitude for sinusoidal strain [] (default=0.2): '))
except ValueError:
    eamp = 0.2
try:
    freq = float(input('frequency for sinusoidal strain [Hz] (default=1.0 Hz): '))
except ValueError:
    freq = 1.0

af = 2*np.pi*freq
s0 = 0                      # [] initial strain

tmax = 5/freq                   # [s] duration time
dt = 1/(40*freq)                # [s] interval time
#dt = 1/(200*freq)                # [s] interval time (for low freq)
t_a = np.arange(0, tmax, dt)    # time after step stress
t_b = np.arange(-2.0,0,dt)      # time before step stress
t = np.concatenate([t_b,t_a])   # whole time 
zeros = np.zeros(len(t_b))
e_a = np.array([eamp*np.sin(af*t) for t in t_a])
e = np.concatenate([zeros,e_a])     # whole strain
e_stat = e[int(0.8*len(e)):]        # 後半部分を抽出（前半は過渡応答を含むから）
# 直列なのでl0=2lとなっている
el = e*2*l                          # [m] elongation

# solution of ODE
sol = odeint(Maxwell_sinuStrain, s0, t_a, args=(eamp,af,E,tau))
s = np.concatenate([zeros,sol[:,0]])        # [MPa] stress
s_stat = s[int(0.8*len(s)):]                # 後半部分を抽出（前半は過渡応答を含むから）
smax = np.max(s_stat)
ind = getNearestIndex2value(s_stat,0)       # 出力信号が0になるindexを抽出
spha = (180/np.pi)*np.arcsin(np.abs(e_stat[ind])/eamp)
integral_s = np.array([s[:k+1].sum()*dt for k in range(len(s))])     # 簡易的なsの積分
e_s = s/E                                   # strain on spring
e_d = integral_s/eta                        # strain on dashpot
el_s = e_s*2*l                              # [m] elongation of spring
el_d = e_d*2*l                              # [m] elongation of dashpot


# scaling for figure
s = s/10**6                     # MPaスケール
smax = smax/10**6               # MPaスケール

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

var_text = r'$\epsilon_{{amp}}$ = {0:.2f}, $f$ = {1:.3f} Hz, $E$ = {2:.1f} MPa, $\eta$ = {3:.1f} kPa s'.format(eamp,freq,E/10**6,eta/10**3)
ax.text(0.4, 0.9, var_text, transform=ax.transAxes)
eq_text = r'd$\sigma$/d$t$ = -$\sigma$/$\tau$ + $E$d$\epsilon$/d$t$'
ax.text(0.4, 0.8, eq_text, transform=ax.transAxes)
res_text = r'$\tau$ = {0:.1f} s, $\omega\tau$ = {1:.3f}'.format(tau, af*tau)
ax.text(0.4, 0.7, res_text, transform=ax.transAxes)
ax.text(0.35, 0.15, '$l_0$', transform=ax.transAxes)
ax.text(0.75, 0.28, '$\epsilon$ (input)', transform=ax.transAxes)
ax.text(0.1, 0.8, '$\omega$ = {0:.2f} s$^{{-1}}$'.format(af), transform=ax.transAxes)
samp_template = r'$\sigma_{{amp}}$ = {0:.2f} MPa'.format(smax)
ax.text(0.75, 0.52, samp_template, transform=ax.transAxes)
spha_template = r'$\theta$ = {0:.1f} $\degree$'.format(spha)
ax.text(0.75, 0.45, spha_template, transform=ax.transAxes)

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

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return rod, point, triangle, rod_da, damper, strain, arrow_p, arrow_n, time_text

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
    # for time
    time_text.set_text(time_template % (t[i]))
    return rod, point, triangle, rod_da, damper, strain, arrow_p, arrow_n, time_text

'''
y_triの中の重要部分は
x_tri1 = np.linspace(a, b,100)
のとき
(xtri - a)/(b - a)
になる 
'''

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f, 
                    init_func=init, blit=True, interval=frame_int, repeat=False)

savefile = "./gif/Maxwell_sinuStrain_ani_(f={0:.3f}Hz).gif".format(freq)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()