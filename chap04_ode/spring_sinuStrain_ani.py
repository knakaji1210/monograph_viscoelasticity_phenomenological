# ordinary differential equation of spring (sinusoidal strain) with animation

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

# external sinusoidal strain
try:
    eamp = float(input('amplitude for sinusoidal strain [] (default=0.25): '))
except ValueError:
    eamp = 0.25
try:
    freq = float(input('frequency for sinusoidal strain [Hz] (default=0.5 Hz): '))
except ValueError:
    freq = 0.5

af = 2*np.pi*freq

# ODE解析で用いる関数の定義
def spring_sinuStrain(s, t, eamp, af, E):
# e: strain, s: stress, E: modulus
# ここではeampとafを指定し、この中でeの関数を作り振動歪みを実現
    e = eamp*np.sin(af*t)           # 振動歪みの関数
    dsdt = E*eamp*af*np.cos(af*t)   # ds/dt = E*de/dt = E*eamp*af*cos(af*t)
    return dsdt

# 1. データ準備
start_time = -2.0   # 開始時間
end_time = 4/freq   # 終了時間
event_time = 0.0    # 振動歪みを加える時刻
time_duration = end_time - start_time  # [s]
time_duration_pre = event_time - start_time
time_duration_post = end_time - event_time
fps = 60    # fps = 30だと足りないので60に変更 
steps = int(time_duration * fps) + 1
interval_ms = 1000 / fps  # 1コマあたりのミリ秒
t = np.linspace(start_time, end_time, steps)
t_pre = t[t < event_time]
t_post = t[t >= event_time]

strain_post = eamp*np.sin(af*t)  # 振動歪みの関数
strain = np.where(t - event_time >= 0, strain_post, 0)

# solution of ODE
s0 = 0                              # 初期条件として定義
sol = odeint(spring_sinuStrain, s0, t_post, args=(eamp,af,E))
stress_pre = np.zeros_like(t_pre)   # ステップ前の応力はゼロ
stress_post = sol[:, 0]             # 応力履歴
stress = np.concatenate([stress_pre, stress_post])

# scaling for figure
e = strain/1.0     # 描画のためのスケーリング
s = stress/10**6   # 描画のためのスケーリング ([MPa]単位に変換)
el = e*l           # [m] elongation

fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111)
ax.grid()
title_text = "spring (Hooke's elasticity): sinusoidal strain"
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

var_text = r'$\epsilon_{{amp}}$ = {0:.2f}, $f$ = {1:.3f} Hz, $E$ = {2:.1f} MPa'.format(eamp,freq,E/10**6)
ax.text(0.4, 0.9, var_text, transform=ax.transAxes)
eq_text = r'$\sigma$ = $E\epsilon$'
ax.text(0.4, 0.8, eq_text, transform=ax.transAxes)
ax.text(0.3, 0.25, '$l_0$', transform=ax.transAxes)
ax.text(0.6, 0.38, '$\sigma$ (output)', transform=ax.transAxes)

# for spring
rod, = ax.plot([],[], 'b', animated=True)
triangle, = ax.plot([],[], 'b', animated=True)
point, = ax.plot([],[], 'ro', markersize='10', animated=True)
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
    return rod, triangle, point, stress, arrow_p, arrow_n, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    x_rod = [3*l/4 + el[i], l + el[i]]
    rod.set_data(x_rod,y_0)
    x_tri = np.linspace(l/4, 3*l/4 + el[i],100)
    y_tri = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri - l/4)/(el[i]+l/2)-np.pi/2+0.1))-1)
    triangle.set_data(x_tri,y_tri)
    point.set_data([l + el[i]],[0])
    a = 0.2    # 見かけ上の振幅
    x_stress = [l, l + a*s[i]/eamp]
    stress.set_data(x_stress,[-1,-1])
    if s[i] > 0:
        arrow_p.set_data([l + a*s[i]/eamp],[-1])
    else:
        arrow_n.set_data([l + a*s[i]/eamp],[-1])
    time_text.set_text(time_template % (t[i]))
    return rod, triangle, point, stress, arrow_p, arrow_n, time_text

# アニメーション実行   
ani = animation.FuncAnimation(fig, update, frames=steps, 
                    init_func=init, blit=True, interval=interval_ms, repeat=False)

savefile = './mp4/spring_sinuStrain_ani_(f={0:.2f}Hz).mp4'.format(freq)
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '60'])

plt.tight_layout()
plt.show()