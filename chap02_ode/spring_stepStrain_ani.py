# ordinary differential equation of spring (step strain) with animation

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib.animation as animation

'''
バネ要素単独では常微分方程式は不要だが、形式的にODEの形で表現してみる
'''

# variables
try:
    E = float(input('modulus [MPa] (default = 0.2 MPa): '))*10**6
except ValueError:
    E = 2*10**5             # [Pa] modulus
l = 0.1                     # [m] equilibrium length
w = 0.5                     # ratio of spring width

# initial condition
try:
    strain_i = float(input('step strain (default = 0.5): '))
except ValueError:
    strain_i = 0.5         # step strain

# ODE解析で用いる関数の定義
def spring_stepStrain(s, t, e, eta):
# e: strain, s: stress, eta: viscosity
# ここでは下でargsとしてs=stress_iを入れてステップ応力を実現
    dsdt = 0    # バネ要素単独では応力は時間変化しないため
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
sol = odeint(spring_stepStrain, s0, t_post, args=(e0,E))
stress_pre = np.zeros_like(t_pre)  # ステップ前の応力はゼロ
stress_post = sol[:, 0]  # 応力履歴
stress = np.concatenate([stress_pre, stress_post])

# scaling for figure
e = strain/1.0     # 描画のためのスケーリング
s = stress/10**6   # 描画のためのスケーリング ([MPa]単位に変換)
el = e*l           # [m] elongation

# 2. グラフの初期設定
fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111)
ax.grid()
title_text = "spring (Hooke's elasticity): step strain"
ax.set_title(title_text)
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
ax.set_xlim(-0.05,0.3)
ax.set_ylim(-5,5)

# for common
y_0 = [0, 0]
ax.plot([0, l/4],y_0, c='b')
ax.plot(0,0, 'ro', markersize='10')
ax.plot([0,l],[-2,-2], c='g')
ax.plot([0,0],[-1.8,-2.2], c='g')
ax.plot([l,l],[-1.8,-2.2], c='g')

var_text = r'$\epsilon_0$ = {0:.2f} MPa, $E$ = {1:.1f} MPa'.format(strain_i,E/10**6)
ax.text(0.5, 0.9, var_text, transform=ax.transAxes)
eq_text = r'$\sigma = E\epsilon_0$'
ax.text(0.5, 0.8, eq_text, transform=ax.transAxes)
ax.text(0.3, 0.25, '$l_0$', transform=ax.transAxes)

# for spring
rod, = ax.plot([],[], 'b', animated=True)
triangle, = ax.plot([],[], 'b', animated=True)
point, = ax.plot([],[], 'ro', markersize='10', animated=True)
# ここでは[],[]としているが、下で***.set_dataで実際の値を入れている

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():
    time_text.set_text('')
    return rod, triangle, point, time_text

# 3. アニメーション更新関数
def update(i):          
    x_rod = [3*l/4 + el[i], l + el[i]]
    rod.set_data(x_rod,y_0)
    x_tri = np.linspace(l/4, 3*l/4 + el[i],100)
    y_tri = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri - l/4)/(el[i]+l/2)-np.pi/2+0.1))-1)
    triangle.set_data(x_tri,y_tri)
    point.set_data([l + el[i]],[0])
    time_text.set_text(time_template % (t[i]))
    return rod, triangle, point, time_text

ani = animation.FuncAnimation(fig, update, frames=steps, 
                    init_func=init, blit=True, interval=interval_ms, repeat=False)

savefile = './mp4/spring_stepStrain_ani.mp4'
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '30'])

plt.show()