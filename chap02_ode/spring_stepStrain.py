# ordinary differential equation of spring (step strain)
'''
バネ要素単独では常微分方程式は不要
'''

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation

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

# 瞬間弾性率の定義
def instantaneous_modulus(t_elapsed, E_i):
    return np.where(t_elapsed >= 0, E_i, 0)

# 1. データ準備
start_time = -2.0   # 開始時間
end_time = 8.0      # 終了時間
time_duration = end_time - start_time  # [s]
fps = 30
steps = int(time_duration * fps) + 1
interval_ms = 1000 / fps  # 1コマあたりのミリ秒
t = np.linspace(start_time, end_time, steps)
t0 = 0.0            # ステップ歪みを加える時刻

# solution of ODE（ここではそれをする必要はない）
stress = strain_i * instantaneous_modulus(t - t0, E)

# scaling for figure
s = stress/10**6                     # MPaスケール

# 2. グラフの初期設定
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

# --- 上段：歪み (Input) ---
ax1.set_xlim(start_time, end_time)
ax1.set_ylim(-0.05, 0.15) # 縦軸を固定
ax1.set_ylabel('Applied Strain, $\epsilon$ /')
ax1.set_title("spring (Hooke's elasticity): step strain")
ax1.grid(True, ls='--')

line_strain, = ax1.plot([], [], color='blue', lw=2, label='Step strain (t = {0:.1f} s)'.format(t0))
ax1.legend(loc='upper right')

# --- 下段：応力 (Response) ---
ax2.set_xlim(start_time, end_time)
ax2.set_ylim(-0.02, 0.04) # 縦軸を固定
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
    curr_t = t[:i]
    # 歪みデータの更新
    strain_history = np.where(curr_t < t0, 0, np.where(curr_t < t0, 0, strain_i))
    line_strain.set_data(curr_t, strain_history)
    
    # 応力データの更新
    line_stress.set_data(t[:i], s[:i])
    return line_strain, line_stress

# アニメーション実行   
ani = animation.FuncAnimation(fig, animate, frames=steps, interval=interval_ms, blit=True)

savefile = './mp4/spring_stepStrain.mp4'
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '30'])

plt.tight_layout()
plt.show()