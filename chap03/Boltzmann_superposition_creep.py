# Fig. 3.6

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# クリープコンプライアンスの定義 (J(t): 単位応力あたりの歪み) ただし粘性流動項なし
def creep_compliance(t_elapsed, J0, J_inf, tau):
    return np.where(t_elapsed >= 0, J0 + (J_inf - J0) * (1 - np.exp(-t_elapsed/tau)), 0)

# 1. データ準備
stress_i = 1.0  # ステップ応力の大きさ
start_time = -2.0  # 開始時間
end_time = 8.0    # 終了時間
time_duration = end_time - start_time  # [s]
fps = 30
steps = int(time_duration * fps) + 1
interval_ms = 1000 / fps  # 1コマあたりのミリ秒
t = np.linspace(start_time, end_time, steps)
t1, t2 = 0.0, 3.0   # ステップ応力を加える時刻
J0 = 2.0    # 瞬間コンプライアンス
J_inf = 4.0 # 定常コンプライアンス
tau = 2.0   # 遅延時間
strain1 = stress_i * creep_compliance(t - t1, J0, J_inf, tau)
strain2 = stress_i * creep_compliance(t - t2, J0, J_inf, tau)
strain_total = strain1 + strain2

stress = np.where(t < t1, 0, np.where(t < t2, stress_i, stress_i * 2))

# 2. グラフの初期設定
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

# --- 上段：応力 (Input) ---
ax1.set_xlim(start_time, end_time)
ax1.set_ylim(-0.2, 2.5) # 縦軸固定
ax1.set_ylabel('Applied Stress, $\sigma$ /MPa')
ax1.set_title('Boltzmann Superposition (Creep) ($J_0$ = {0:.1f} MPa$^{{-1}}$, $J_\\infty$ = {1:.1f} MPa$^{{-1}}$, $\\tau$ = {2:.1f} s)'.format(J0, J_inf, tau))
ax1.grid(True, ls='--')
line_stress, = ax1.step([], [], where='post', color='red', lw=2, label='Two-step stress (t = {0:.1f} s)'.format(t1))
ax1.legend(loc='upper right')

# --- 下段：歪み (Response) ---
ax2.set_xlim(start_time, end_time)
ax2.set_ylim(-0.5, 10.0) # 縦軸固定
ax2.set_xlabel('Time /s')
ax2.set_ylabel('Strain, $\epsilon$ /')
ax2.grid(True, ls='--')

line_total, = ax2.plot([], [], 'b-', lw=2.5, label='Total Strain (Sum)')
line_str1, = ax2.plot([], [], 'g--', alpha=0.5, label='Response to Step 1 (t = {0:.1f} s)'.format(t1))
line_str2, = ax2.plot([], [], 'y--', alpha=0.5, label='Response to Step 1 (t = {0:.1f} s)'.format(t2))
ax2.legend(loc='upper left')

# 3. アニメーション更新関数
def animate(i):
    # 応力データの更新
    line_stress.set_data(t[:i], stress[:i])
    
    # 歪みデータの更新
    line_str1.set_data(t[:i], strain1[:i])
    line_str2.set_data(t[:i], strain2[:i])
    line_total.set_data(t[:i], strain_total[:i])
    return line_stress, line_total, line_str1, line_str2

# アニメーション実行
ani = animation.FuncAnimation(fig, animate, frames=steps, interval=interval_ms, blit=True)

savefile = './mp4/Boltzmann_creep.mp4'
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '30'])

plt.tight_layout()
plt.show()
