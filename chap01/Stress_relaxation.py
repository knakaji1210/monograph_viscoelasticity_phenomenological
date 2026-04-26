# Figure 1.8

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# 緩和弾性率の定義
def relaxation_modulus(t_elapsed, G_inf, G_i, tau):
    return np.where(t_elapsed >= 0, G_inf + (G_i - G_inf) * np.exp(-t_elapsed/tau), 0)

# 1. データ準備
strain_i = 1.0  # ステップ歪みの大きさ
t = np.linspace(-2, 10, 400)
t0 = 0.0        # ステップ歪みを加える時刻
G_inf = 1.0     # 平衡弾性率
G_i = 3.0       # 瞬間弾性率
tau = 3.0       # 緩和時間
stress = strain_i * relaxation_modulus(t - t0, G_inf, G_i, tau)

# 2. グラフの初期設定
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

# --- 上段：歪み (Input) ---
ax1.set_xlim(-2, 10)
ax1.set_ylim(-0.2, 2.0) # 縦軸を固定
ax1.set_ylabel('Applied Strain, $\epsilon$ /')
ax1.set_title('Strees Relaxation ($G_i$ = {0:.1f} MPa, $G_\\infty$ = {1:.1f} MPa, $\\tau$ = {2:.1f} s)'.format(G_i, G_inf, tau))
ax1.grid(True, ls='--')

line_strain, = ax1.step([], [], where='post', color='blue', lw=2, label='Step strain (t = {0:.1f} s)'.format(t0))
ax1.legend(loc='upper right')

# --- 下段：応力 (Response) ---
ax2.set_xlim(-2, 10)
ax2.set_ylim(-0.4, 4.0) # 縦軸を固定
ax2.set_xlabel('Time /s')
ax2.set_ylabel('Stress, $\sigma$ /MPa')
ax2.grid(True, ls='--')

line_stress, = ax2.plot([], [], color='red', lw=2, label='Response to step strain (t = {0:.1f} s)'.format(t0))
ax2.legend(loc='upper right')

# 3. アニメーション更新関数
def animate(i):
    curr_t = t[:i]
    # 歪みデータの更新
    strain_history = np.where(curr_t < t0, 0, np.where(curr_t < t0, 0, strain_i))
    line_strain.set_data(curr_t, strain_history)
    
    # 応力データの更新
    line_stress.set_data(t[:i], stress[:i])
    return line_strain, line_stress

# アニメーション実行
ani = animation.FuncAnimation(fig, animate, frames=len(t), interval=25, blit=True)

ani.save('./gif/Stress_relaxation.gif', dpi=300)

plt.tight_layout()
plt.show()
