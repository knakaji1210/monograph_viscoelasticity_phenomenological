# Fig. 3.2

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# 緩和弾性率の定義
def relaxation_modulus(t_elapsed, E_inf, E_i, tau):
    return np.where(t_elapsed >= 0, E_inf + (E_i - E_inf) * np.exp(-t_elapsed/tau), 0)

# 1. データ準備
strain_i = 1.0  # ステップ歪みの大きさ
start_time = -2.0  # 開始時間
end_time = 8.0    # 終了時間
time_duration = end_time - start_time  # [s]
fps = 30
steps = int(time_duration * fps) + 1
interval_ms = 1000 / fps  # 1コマあたりのミリ秒
t = np.linspace(start_time, end_time, steps)
t1, t2 = 0.0, 2.0   # ステップ歪みを加える時刻
E_inf = 0.5         # 平衡弾性率
E_i = 3.0           # 瞬間弾性率
tau = 2.0           # 緩和時間
stress1 = strain_i * relaxation_modulus(t - t1, E_inf, E_i, tau)
stress2 = strain_i * relaxation_modulus(t - t2, E_inf, E_i, tau)
stress_total = stress1 + stress2

strain = np.where(t < t1, 0, np.where(t < t2, 1.0, 2.0))

# 2. グラフの初期設定
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

# --- 上段：歪み (Input) ---
ax1.set_xlim(start_time, end_time)
ax1.set_ylim(-0.2, 3) # 縦軸を固定
ax1.set_ylabel('Applied Strain, $\epsilon$ /')
ax1.set_title('Boltzmann Superposition (Strees Relaxation) ($E_i$ = {0:.1f} MPa, $E_\\infty$ = {1:.1f} MPa, $\\tau$ = {2:.1f} s)'.format(E_i, E_inf, tau))

ax1.grid(True, ls='--')
line_strain, = ax1.step([], [], where='post', color='blue', lw=2)

# --- 下段：応力 (Response) ---
ax2.set_xlim(start_time, end_time)
ax2.set_ylim(-0.2, 5.0) # 縦軸を固定
ax2.set_xlabel('Time /s')
ax2.set_ylabel('Stress, $\sigma$ /MPa')
ax2.grid(True, ls='--')

line_s1, = ax2.plot([], [], 'y--', alpha=0.5, label='Response to Step 1 (t = {0:.1f} s)'.format(t1))
line_s2, = ax2.plot([], [], 'g--', alpha=0.5, label='Response to Step 2 (t = {0:.1f} s)'.format(t2))
line_total, = ax2.plot([], [], 'r-', lw=2.5, label='Total Stress (Sum)')
ax2.legend(loc='upper right')

# 3. アニメーション更新関数
def animate(i):
     # 歪みデータの更新
    line_strain.set_data(t[:i], strain[:i])
    
    # 応力データの更新
    line_s1.set_data(t[:i], stress1[:i])
    line_s2.set_data(t[:i], stress2[:i])
    line_total.set_data(t[:i], stress_total[:i])
    return line_strain, line_s1, line_s2, line_total

# アニメーション実行
ani = animation.FuncAnimation(fig, animate, frames=steps, interval=interval_ms, blit=True)

savefile = './mp4/Boltzmann_relax_stepStrain.mp4'
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '30'])

plt.tight_layout()
plt.show()
