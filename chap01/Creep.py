# Figure 1.10

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# クリープコンプライアンスの定義 (J(t): 単位応力あたりの歪み)
def creep_compliance(t_elapsed, J0, J_inf, eta, tau):
    return np.where(t_elapsed >= 0, J0 + eta * t_elapsed + (J_inf - J0) * (1 - np.exp(-t_elapsed/tau)), 0)

# 1. データ準備
stress_i = 1.0  # ステップ応力の大きさ
t = np.linspace(-2, 10, 500)
t0 = 0.0    # 荷重（応力）を加える時刻
J0 = 2.0    # 瞬間コンプライアンス
J_inf = 4.0 # 定常コンプライアンス
eta = 1.0   # 粘性係数（本来はtau/J0である）
tau = 2.0   # 遅延時間

# 各荷重による個別の歪み応答
strain1 = stress_i * creep_compliance(t - t0, J0, J_inf, eta, tau)
strain2 = stress_i * creep_compliance(t - t0, J0, J_inf, 0, tau)    # 粘性流動項なしの応答
strain3 = stress_i * creep_compliance(t - t0, J0, J0, eta, tau)     # 緩和項なしの応答 (粘性流動のみ)

# 2. グラフの初期設定
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

# --- 上段：応力 (Input) ---
ax1.set_xlim(-2, 10)
ax1.set_ylim(-0.2, 2.0) # 縦軸固定
ax1.set_ylabel('Applied Stress, $\sigma$ /MPa')
ax1.set_title('Creep ($J_0$ = {0:.1f} MPa$^{{-1}}$, $J_\\infty$ = {1:.1f} MPa$^{{-1}}$, $\\eta$ = {2:.1f} MPa$\cdot$s, $\\tau$ = {3:.1f} s)'.format(J0, J_inf, eta, tau))
ax1.grid(True, ls='--')
line_stress, = ax1.step([], [], where='post', color='blue', lw=2, label='Step stress (t = {0:.1f} s)'.format(t0))
ax1.legend(loc='upper right')

# --- 下段：歪み (Response) ---
ax2.set_xlim(-2, 10)
ax2.set_ylim(-0.2, 15.0) # 縦軸固定
ax2.set_xlabel('Time /s')
ax2.set_ylabel('Strain, $\epsilon$ /')
ax2.grid(True, ls='--')

line_stress1, = ax2.plot([], [], 'r-', lw=2, label='total creep response (with viscous flow)')
line_stress2, = ax2.plot([], [], 'g--', lw=2, label='creep response (no viscous flow)')
line_stress3, = ax2.plot([], [], 'b--', lw=2, label='creep response (no relaxation term)')
ax2.legend(loc='upper left')

# 3. アニメーション更新関数
def animate(i):
    curr_t = t[:i]
    # 応力履歴の更新
    stress_history = np.where(curr_t < t0, 0, np.where(curr_t < t0, 0, stress_i))
    line_stress.set_data(curr_t, stress_history)
    
    # 歪みデータの更新
    line_stress1.set_data(t[:i], strain1[:i])
    line_stress2.set_data(t[:i], strain2[:i])
    line_stress3.set_data(t[:i], strain3[:i])
    return line_stress, line_stress1, line_stress2, line_stress3

# アニメーション実行
ani = animation.FuncAnimation(fig, animate, frames=len(t), interval=20, blit=True)

ani.save('./gif/Creep.gif', dpi=300)

plt.tight_layout()
plt.show()
