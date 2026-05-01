# ordinary differential equation of dashpot (step stress)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation

'''
テキストの式(2.6)をベースに組み立てる
'''

# variables
try:
    eta = float(input('viscosity [kPa s] (default = 500.0 kPa s): '))*10**3
except ValueError:
    eta = 5*10**5             # [Pa s] viscosity

# initial condition
try:
    stress_i = float(input('step stress [MPa] (default = 0.02 MPa): '))*10**6
except ValueError:
    stress_i = 0.02*10**6         # [Pa] step stress

# ODE解析で用いる関数の定義
def dashpot_stepStress(e, t, s, eta):
# e: strain, s: stress, eta: viscosity
# ここでは下でargsとしてs0=stress_iを入れてステップ応力を実現
    dedt = s/eta    # (2.6)
    return dedt

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

stress = np.where(t - event_time >= 0, stress_i, 0)

# solution of ODE
s0 = stress_i   # ODEの引数として入れるためにこの形で定義
e0 = 0.0        # ステップ応力を加える前の歪みはゼロとするため、初期条件として定義
sol = odeint(dashpot_stepStress, e0, t_post, args=(s0,eta))
strain_pre = np.zeros_like(t_pre)  # ステップ前の歪みはゼロ
strain_post = sol[:, 0]     # 歪み履歴
strain = np.concatenate([strain_pre, strain_post])

# scaling for figure
e = strain/1.0     # 描画のためのスケーリング
s = stress/10**6   # 描画のためのスケーリング ([MPa]単位に変換)

# 2. グラフの初期設定
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

# --- 上段：応力 (Input) ---
ax1.set_xlim(start_time, end_time)
ax1.set_ylim(-np.max(s)*0.5, np.max(s)*1.5) # 縦軸を固定
ax1.set_ylabel('Applied stress, $\sigma$ /MPa')
ax1.set_title("dashpot (Newton's viscosity): step stress")
ax1.grid(True, ls='--')

line_stress, = ax1.plot([], [], color='red', lw=2, label='Step stress (t = {0:.1f} s)'.format(event_time))
ax1.legend(loc='upper right')

# --- 下段：歪み (Response) ---
ax2.set_xlim(start_time, end_time)
ax2.set_ylim(-np.max(e)*0.1, np.max(e)*1.5) # 縦軸を固定
ax2.set_xlabel('$t$ /s')
ax2.set_ylabel('Strain, $\epsilon$ /')
ax2.grid(True, ls='--')

line_strain, = ax2.plot([], [], color='blue', lw=2, label='Response to step stress')
ax2.legend(loc='upper right')

var_text = r'$\sigma_0$ = {0:.2f} MPa, $\eta$ = {1:.1f} kPa s'.format(stress_i/10**6,eta/10**3)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'd$\epsilon$/d$t$ = $\sigma_0$/$\eta$'
ax2.text(0.1, 0.9, eq_text, transform=ax2.transAxes)

# 3. アニメーション更新関数
def animate(i):
    # 応力データの更新
    line_stress.set_data(t[:i], s[:i])
    
    # 歪みデータの更新
    line_strain.set_data(t[:i], e[:i])
    return line_stress, line_strain

ani = animation.FuncAnimation(fig, animate, frames=steps, interval=interval_ms, blit=True)

savefile = './mp4/dashpot_stepStress.mp4'
#ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '30'])

plt.tight_layout()
plt.show()