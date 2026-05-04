# Figure 1.8

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# 緩和弾性率の定義
def relaxation_modulus(t_elapsed, E_inf, E_i, tau):
    return np.where(t_elapsed >= 0, E_inf + (E_i - E_inf) * np.exp(-t_elapsed/tau), 0)

# データ準備
strain_i = 1.0      # ステップ歪みの大きさ
start_time = -2.0   # [s] 開始時間
end_time = 8.0      # [s] 終了時間
time_duration = end_time - start_time   # [s]　継続時間
fps = 30            # 1秒あたりのフレーム数 
steps = int(time_duration * fps) + 1    # 総フレーム数
interval_ms = 1000 / fps                # 1コマあたりのミリ秒
t = np.linspace(start_time, end_time, steps)
t0 = 0.0        # ステップ歪みを加える時刻
E_inf = 1.0     # 平衡弾性率
E_i = 3.0       # 瞬間弾性率
tau = 2.5       # 緩和時間
stress1 = strain_i * relaxation_modulus(t - t0, E_inf, E_i, tau)
stress2 = strain_i * relaxation_modulus(t - t0, E_inf, E_inf, tau)

# グラフの初期設定
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

# --- 上段：歪み (Input) ---
ax1.set_xlim(start_time, end_time)
ax1.set_ylim(-0.2, 2.0) # 縦軸を固定
ax1.set_ylabel('Applied Strain, $\epsilon$ /')
ax1.set_title('Strees Relaxation ($E_i$ = {0:.1f} MPa, $E_\\infty$ = {1:.1f} MPa, $\\tau$ = {2:.1f} s)'.format(E_i, E_inf, tau))
ax1.grid(True, ls='--')

line_strain, = ax1.step([], [], where='post', color='blue', lw=2, label='Step strain (t = {0:.1f} s)'.format(t0))
ax1.legend(loc='upper right')

# --- 下段：応力 (Response) ---
ax2.set_xlim(start_time, end_time)
ax2.set_ylim(-0.4, 4.0) # 縦軸を固定
ax2.set_xlabel('Time /s')
ax2.set_ylabel('Stress, $\sigma$ /MPa')
ax2.grid(True, ls='--')

line_stress1, = ax2.plot([], [], color='red', lw=2, label='Response to step strain (total response)')
line_stress2, = ax2.plot([], [], '--', color='green', lw=2, label='Response to step strain ($E_\\infty$ only)')
ax2.legend(loc='upper right')

# アニメーション更新関数
def animate(i):
    curr_t = t[:i]
    # 歪みデータの更新
    strain_history = np.where(curr_t < t0, 0, np.where(curr_t < t0, 0, strain_i))
    line_strain.set_data(curr_t, strain_history)
    
    # 応力データの更新
    line_stress1.set_data(t[:i], stress1[:i])
    line_stress2.set_data(t[:i], stress2[:i])
    return line_strain, line_stress1, line_stress2

# アニメーション実行   
ani = animation.FuncAnimation(fig, animate, frames=steps, interval=interval_ms, blit=True)

savefile = './mp4/Stress_relaxation.mp4'
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '30'])

plt.tight_layout()
plt.show()
