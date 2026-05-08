# ordinary differential equation of Voigt model (step stress) with animation

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.animation import FuncAnimation

'''
テキストの式(1.22)をベースに組み立てる
'''

def Voigt_sinuStress(e, t, samp, af, tau):
# e: strain, s: stress, tau: retardation time
# ここではsampとafを指定し、この中でsの関数を作り振動応力を実現
    s = samp*np.sin(af*t)
    dedt = (s/E - e)/tau        # (1.21)
    return dedt

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
tau = eta/E                     # [s] retardation time
l = 0.1                     # [m] equilibrium length
w = 0.5                     # ratio of dashpot width

# external sinusoidal strain
try:
    samp = float(input('amplitude for sinusoidal stress [MPa] (default=0.05): '))*10**6
except ValueError:
    samp = 0.05*10**6
try:
    T = float(input('peirod for sinusoidal strain [s] (default=2.0): '))
except ValueError:
    T = 2.0

af = 2*np.pi/T
e0 = 0                      # [] initial strain

tmax = 5*T                  # [s] duration time
dt = 0.1                    # [s] interval time
t_a = np.arange(0, tmax, dt)    # time after step stress
t_b = np.arange(-2.0,0,dt) # time before step stress
t = np.concatenate([t_b,t_a])   # whole time 
zeros = np.zeros(len(t_b))
s_a = np.array([samp*np.sin(af*t) for t in t_a])
s = np.concatenate([zeros,s_a])     # whole stress
s_stat = s[int(0.8*len(s)):]        # 後半部分を抽出（前半は過渡応答を含むから）

# solution of ODE
sol = odeint(Voigt_sinuStress, e0, t_a, args=(samp,af,tau))
e = np.concatenate([zeros,sol[:,0]])        # [] strain
e_stat = e[int(0.8*len(e)):]                # 後半部分を抽出（前半は過渡応答を含むから）
emax = np.max(e_stat)
ind = getNearestIndex2value(e_stat,0)       # 出力信号が0になるindexを抽出
epha = (180/np.pi)*np.arcsin(np.abs(s_stat[ind])/samp) 
dedt = np.array([0.0]+[(e[k+1]-e[k])/(t[k+1]-t[k]) for k in range(len(e)-1)])     # 簡易的なeの微分
s_s = E*e                                   # stress on spring
s_d = eta*dedt                              # stress on dashpot
# 並列なのでl0=lとなっている
el = e*l                                    # [m] elongation

# scaling for figure
s = s/10**6                     # MPaスケール
s_s = s_s/10**6                 # MPaスケール
s_d = s_d/10**6                 # MPaスケール

fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111, xlabel='$t$ /s')
ax.grid()
title_text = "Voigt model: sinusoidal stress"
ax.set_title(title_text)
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
ax.set_xlim(-0.05,0.3)
ax.set_ylim(-5,5)
# for common
y_0 = [0, 0]
y_1 = [1, 1]
y_2 = [-1, -1]
ax.plot([0, l/4],y_0, c='b')
ax.plot([l/4, 5*l/4],[-3,-3], c='g')
ax.plot([l/4, l/4],[-2.8,-3.2], c='g')
ax.plot([5*l/4, 5*l/4],[-2.8,-3.2], c='g')
ax.plot([l/4,l/4],[1,-1], c='b')
ax.plot([l/4, (0.2+1/4)*l],y_1, c='b')
ax.plot([l/4, l/2],y_2, c='b')
ax.plot(0,0,'ro', markersize='10')
# for dashpot
x_d1 = [(0.2+1/4)*l, (0.8+1/4)*l]
y_d1 = [w+1, w+1]
y_d2 = [-w+1, -w+1]
ax.plot(x_d1,y_d1, c='b')
ax.plot(x_d1,y_d2, c='b')
ax.plot([(0.2+1/4)*l,(0.2+1/4)*l],[w+1,-w+1], c='b')
rect = patches.Rectangle(xy=((0.2+1/4)*l, -w+1), width=0.6*l, height=2*w, facecolor='y')
ax.add_patch(rect)

var_text = r'$\sigma_{{amp}}$ = {0:.2f} MPa, $T$ = {1:.1f} s, $E$ = {2:.1f} MPa, $\eta$ = {3:.1f} kPa s'.format(samp/10**6,T,E/10**6,eta/10**3)
ax.text(0.4, 0.9, var_text, transform=ax.transAxes)
eq_text = r'd$\epsilon$/d$t$ = ($\sigma$/$E$ - $\epsilon$)/$\tau$'
ax.text(0.4, 0.8, eq_text, transform=ax.transAxes)
res_text = r'$\tau$ = {0:.1f} s'.format(tau)
ax.text(0.4, 0.7, res_text, transform=ax.transAxes)
ax.text(0.35, 0.15, '$l_0$', transform=ax.transAxes)
ax.text(0.75, 0.28, '$\sigma$ (input)', transform=ax.transAxes)
ax.text(0.1, 0.8, '$\omega_f$ = {0:.2f} s$^{{-1}}$'.format(af), transform=ax.transAxes)
eamp_template = r'$\epsilon_{{amp}}$ = {0:.2f}'.format(emax)
ax.text(0.75, 0.52, eamp_template, transform=ax.transAxes)
epha_template = r'$\theta$ = {0:.1f} $\degree$'.format(epha)
ax.text(0.75, 0.45, epha_template, transform=ax.transAxes)

# for common
bar, = ax.plot([],[], 'b', animated=True)
rod, = ax.plot([],[], 'b', animated=True)
point, = ax.plot([],[], 'ro', markersize='10', animated=True)
# for spring
rod_sp, = ax.plot([],[], 'b', animated=True)
triangle, = ax.plot([],[], 'b', animated=True)
# for dashpot
rod_da, = ax.plot([],[], 'b', animated=True)
damper, = ax.plot([],[], 'b', lw=4, animated=True)

# ここでは[],[]としているが、下で***.set_dataで実際の値を入れている

# for stress
stress, = ax.plot([],[], 'r', lw=2, animated=True)
arrow_p, = ax.plot([],[], 'r', marker=9, markersize='10', animated=True)
arrow_n, = ax.plot([],[], 'r', marker=8, markersize='10', animated=True)

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return bar, rod, point, rod_sp, triangle, rod_da, damper, stress, arrow_p, arrow_n, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    # for common
    x_o = l + el[i]
    x_c = l/4 + x_o
    bar.set_data([x_c,x_c],[1,-1])
    x_rod = [x_c, l/4+x_c]
    rod.set_data(x_rod,y_0)
    point.set_data([l/4+x_c],[0])
    # for spring
    x_rod_sp = [x_o, x_c]
    rod_sp.set_data(x_rod_sp,y_2)
    x_tri = np.linspace(l/2, x_o,100)
    y_tri = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri - l/2)/(el[i]+l/2)-np.pi/2+0.1))-3)
    triangle.set_data(x_tri,y_tri)
    # for dashpot
    x_rod_da = [x_o-l/4, x_c]
    x_damp = x_o-l/4
    y_damp = 0.7*w
    x_damper = [x_damp, x_damp]
    y_damper = [y_damp+1, -y_damp+1]
    rod_da.set_data(x_rod_da,y_1)
    damper.set_data(x_damper,y_damper)
    a = 0.2    # 見かけ上の振幅
    x_stress = [3*l/2, 3*l/2 + a*s[i]]
    stress.set_data(x_stress,[-2,-2])
    if s[i] > 0:
        arrow_p.set_data([3*l/2 + a*s[i]],[-2])
    else:
        arrow_n.set_data([3*l/2 + a*s[i]],[-2])
    time_text.set_text(time_template % (t[i]))
    return bar, rod, point, rod_sp, triangle, rod_da, damper, stress, arrow_p, arrow_n, time_text

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

savefile = "./gif/Voigt_sinuStress_(sigma={0:.2f}MPa,T={1:.1f}s,mod={2:.1f}MPa,eta={3:.1f}kPas).gif".format(samp/10**6,T,E/10**6,eta/10**3)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()