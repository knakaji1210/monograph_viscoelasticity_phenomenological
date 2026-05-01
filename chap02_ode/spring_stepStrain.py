# ordinary differential equation of spring (step strain)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation

'''
バネ要素単独では常微分方程式は不要だが、形式的にODEの形で表現してみる
'''

# variables
try:
    E = float(input('modulus [MPa] (default = 0.2 MPa): '))*10**6
except ValueError:
    E = 2*10**5             # [Pa] modulus

# initial condition
try:
    strain_i = float(input('step strain (default = 0.1): '))
except ValueError:
    strain_i = 0.1         # step strain

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

# 2. グラフの初期設定
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

# --- 上段：歪み (Input) ---
ax1.set_xlim(start_time, end_time)
ax1.set_ylim(-np.max(e)*0.5, np.max(e)*1.5) # 縦軸を固定
ax1.set_ylabel('Applied Strain, $\epsilon$ /')
ax1.set_title("spring (Hooke's elasticity): step strain")
ax1.grid(True, ls='--')

line_strain, = ax1.plot([], [], color='blue', lw=2, label='Step strain (t = {0:.1f} s)'.format(event_time))
ax1.legend(loc='upper right')

# --- 下段：応力 (Response) ---
ax2.set_xlim(start_time, end_time)
ax2.set_ylim(-np.max(s)*0.5, np.max(s)*1.5) # 縦軸を固定
ax2.set_xlabel('$t$ /s')
ax2.set_ylabel('Stress, $\sigma$ /MPa')
ax2.grid(True, ls='--')

line_stress, = ax2.plot([], [], color='red', lw=2, label='Response to step strain')
ax2.legend(loc='upper right')

var_text = r'$\epsilon_0$ = {0:.2f}, $E$ = {1:.1f} MPa'.format(strain_i, E/10**6)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'$\sigma = E\epsilon_0$'
ax2.text(0.1, 0.9, eq_text, transform=ax2.transAxes)

# 3. アニメーション更新関数
def animate(i):
    # 歪みデータの更新
    line_strain.set_data(t[:i], e[:i])

    # 応力データの更新
    line_stress.set_data(t[:i], s[:i])

    return line_strain, line_stress

# アニメーション実行   
ani = animation.FuncAnimation(fig, animate, frames=steps, interval=interval_ms, blit=True)

savefile = './mp4/spring_stepStrain.mp4'
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '30'])

plt.tight_layout()
plt.show()