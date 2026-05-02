# ordinary differential equation of Maxwell model (step strain) with animation

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib.animation as animation

'''
テキストの式(2.12)をベースに組み立てる
'''

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
l = 0.1                         # [m] equilibrium length
w = 0.5                         # ratio of dashpot width

# initial condition
try:
    strain_i = float(input('step strain [] (default = 0.2): '))
except ValueError:
    strain_i = 0.2               # [] step strain

# ODE解析で用いる関数の定義
def Maxwell_stepStrain(s, t, tau):
# e: strain, s: stress, tau: retardation time
# ステップ歪みe0=strain_iを加えるので、式(2.11)のde/dtの項が0となっている
    dsdt = -s/tau               # (2.12)
    return dsdt

# 1. データ準備
start_time = -2.0   # 開始時間
end_time = 8.0      # 終了時間
event_time = 0.0    # ステップ歪みを加える時刻
time_duration = end_time - start_time  # [s]
time_duration_pre = event_time - start_time
time_duration_post = end_time - event_time
fps = 30
steps = int(time_duration * fps) + 1
interval_ms = 1000 / fps  # 1コマあたりのミリ秒
t = np.linspace(start_time, end_time, steps)
t_pre = t[t < event_time]
t_post = t[t >= event_time]

strain = np.where(t - event_time >= 0, strain_i, 0)

# solution of ODE
e0 = strain_i       # ODEの引数として入れるためにこの形で定義
s0 = E * strain_i   # バネ要素の応力はステップ歪みに対して即座に応答するため、初期条件として定義
sol = odeint(Maxwell_stepStrain, s0, t_post, args=(tau,))
stress_pre = np.zeros_like(t_pre)  # ステップ前の応力はゼロ
stress_post = sol[:, 0]  # 応力履歴
stress = np.concatenate([stress_pre, stress_post])

dt = t[1] - t[0]   # 時間刻み
integral_stress = np.array([stress[:k+1].sum()*dt for k in range(len(stress))])     # 簡易的なsの積分

# scaling for figure
e = strain/1.0     # 描画のためのスケーリング
s = stress/10**6   # 描画のためのスケーリング ([MPa]単位に変換)
e_s = s/(E/10**6)  # バネの歪み
e_d = integral_stress/eta   # ダッシュポットの歪み
# 直列なのでl0=2*lとなっている
el = e*2*l         # [m] 全体の伸び
el_s = e_s*2*l     # [m] バネの伸び
el_d = e_d*2*l     # [m] ダッシュポットの伸び

# 2. グラフの初期設定
fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111, xlabel='$t$ /s')
ax.grid()
title_text = "Maxwell model: step strain"
ax.set_title(title_text)
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
ax.set_xlim(-0.05,0.4)
ax.set_ylim(-5,5)
# for common
y_0 = [0, 0]
ax.plot([0, 0.08*l],y_0, c='b')
ax.plot(0,0,'ro', markersize='10')
ax.plot([0,2*l],[-2,-2], c='g')
ax.plot([0,0],[-1.8,-2.2], c='g')
ax.plot([2*l,2*l],[-1.8,-2.2], c='g')
# for dashpot
x_d1 = [0.08*l, 0.92*l]
y_d1 = [w, w]
y_d2 = [-w, -w]
ax.plot(x_d1,y_d1, c='b')
ax.plot(x_d1,y_d2, c='b')
ax.plot([0.08*l,0.08*l],[w,-w], c='b')
rect = patches.Rectangle(xy=(0.08*l, -w), width=0.83*l, height=2*w, facecolor='y')
ax.add_patch(rect)

var_text = r'$\epsilon_0$ = {0:.1f}, $E$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s'.format(strain_i,E/10**6,eta/10**3)
ax.text(0.5, 0.9, var_text, transform=ax.transAxes)
eq_text = r'd$\sigma$/d$t$ = -$\sigma$/$\tau$'
ax.text(0.5, 0.8, eq_text, transform=ax.transAxes)
res_text = r'$\tau$ = {0:.1f} s'.format(tau)
ax.text(0.5, 0.7, res_text, transform=ax.transAxes)
ax.text(0.35, 0.25, '$l_0$', transform=ax.transAxes)

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

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return rod, point, triangle, rod_da, damper, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    # for common
    x_o = 7*l/4 + el[i]
    x_rod = [x_o, l/4+x_o]
    rod.set_data(x_rod,y_0)
    point.set_data([l/4+x_o],[0])
    # for dashpot
    x_da = l/2 + el_d[i]
    x_rod_da = [x_da, x_da+3*l/4]
    x_damp = x_da
    y_damp = 0.7*w
    x_damper = [x_damp, x_damp]
    y_damper = [y_damp, -y_damp]
    rod_da.set_data(x_rod_da,y_0)
    damper.set_data(x_damper,y_damper)
    time_text.set_text(time_template % (t[i]))
    # for spring
    x_sp = x_o-el_s[i]-l/2
    x_tri = np.linspace(x_sp, x_o,100)
    y_tri = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri - x_sp)/(x_o-x_sp)-np.pi/2+0.1))-1)
    triangle.set_data(x_tri,y_tri)
    return rod, point, triangle, rod_da, damper, time_text

'''
y_triの中の重要部分は
x_tri1 = np.linspace(a, b,100)
のとき
(xtri - a)/(b - a)
になる 
'''

ani = animation.FuncAnimation(fig, update, frames=steps, 
                    init_func=init, blit=True, interval=interval_ms, repeat=False)

savefile = "./mp4/Maxwell_stepStrain_ani.mp4"
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '30'])

plt.show()