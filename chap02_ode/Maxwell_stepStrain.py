# ordinary differential equation of Maxwell model (step strain)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
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

# 2. グラフの初期設定
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

# --- 上段：歪み (Input) ---
ax1.set_xlim(start_time, end_time)
#ax1.set_ylim(-np.max(e)*0.5, np.max(e)*1.5) # 縦軸を固定
ax1.set_ylim(-np.max(e_d)*0.1, np.max(e_d)*1.5) # 縦軸を固定
ax1.set_ylabel('Applied Strain, $\epsilon$ /')
ax1.set_title("Maxwell model: step strain")
ax1.grid(True, ls='--')

line_strain, = ax1.plot([], [], color='blue', lw=2, label='Step strain (t = {0:.1f} s)'.format(event_time))
line_strain_s, = ax1.plot([], [], color='green', ls="dashed", lw=1, label='$\epsilon$ (spring)')
line_strain_d, = ax1.plot([], [], color='orange', ls="dashed", lw=1, label='$\epsilon$ (dashpot)')
ax1.legend(loc='upper right')

# --- 下段：応力 (Response) ---
ax2.set_xlim(start_time, end_time)
ax2.set_ylim(-np.max(s)*0.5, np.max(s)*1.5) # 縦軸を固定
ax2.set_xlabel('$t$ /s')
ax2.set_ylabel('Stress, $\sigma$ /MPa')
ax2.grid(True, ls='--')

line_stress, = ax2.plot([], [], color='red', lw=2, label='Response to step strain')
ax2.legend(loc='upper right')

var_text = r'$\epsilon_0$ = {0:.1f}, $E$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s'.format(e0,E/10**6,eta/10**3)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'd$\sigma$/d$t$ = -$\sigma$/$\tau$'
ax2.text(0.1, 0.9, eq_text, transform=ax2.transAxes)
res_text = r'$\tau$ = {0:.1f} s'.format(tau)
ax2.text(0.1, 0.8, res_text, transform=ax2.transAxes)

# 3. アニメーション更新関数
def animate(i):
    # 歪みデータの更新
    line_strain.set_data(t[:i], e[:i])
    line_strain_s.set_data(t[:i], e_s[:i])
    line_strain_d.set_data(t[:i], e_d[:i])

    # 応力データの更新
    line_stress.set_data(t[:i], s[:i])

    return line_strain, line_strain_s, line_strain_d, line_stress

# アニメーション実行   
ani = animation.FuncAnimation(fig, animate, frames=steps, interval=interval_ms, blit=True)

savefile = './mp4/Maxwell_stepStrain.mp4'
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '30'])

plt.tight_layout()
plt.show()