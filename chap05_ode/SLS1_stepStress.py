# SLS Iモデルの常微分方程式（ステップ応力）

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation

'''
テキストの式(5.12)をベースに組み立てる
'''

# 変数の設定
try:
    E1 = float(input('modulus 1 [MPa] (default = 1.0 MPa): '))*10**6
except ValueError:
    E1 = 10**6                  # [Pa] 弾性率
try:
    E2 = float(input('modulus 2 [MPa] (default = 0.2 MPa): '))*10**6
except ValueError:
    E2 = 2*10**5                # [Pa] 弾性率
try:
    eta = float(input('viscosity [kPa s] (default = 500.0 kPa s): '))*10**3
except ValueError:
    eta = 5*10**5               # [Pa s] 粘度

insMod = E1                 # [Pa] 瞬間弾性率
infMod = E1*E2/(E1+E2)      # [Pa] 緩和弾性率
k = insMod/infMod
tau = eta/E2                # [s] 遅延時間

# 初期条件の設定
try:
    stress_i = float(input('step stress [MPa] (default = 0.04 MPa): '))*10**6
except ValueError:
    stress_i = 0.04*10**6         # [Pa] ステップ応力

# ODE解析で用いる関数の定義
def SLS1_stepStress(e, t, s, infMod, tau):
# e: 歪み, s: 応力, infMod: 緩和弾性率, tau: 遅延時間
# ここでは下でargsとしてs=stress_iを入れてステップ応力を実現
    dedt = (s/infMod - e)/tau   # (5.12)
    return dedt

# データ準備
start_time = -2.0   # 開始時間
end_time = 8.0      # 終了時間
event_time = 0.0    # ステップ歪みを加える時刻
time_duration = end_time - start_time       # [s] 継続時間
time_duration_pre = event_time - start_time # [s] ステップ前の継続時間
time_duration_post = end_time - event_time  # [s] ステップ後の継続時間
fps = 30            # 1秒あたりのフレーム数
steps = int(time_duration * fps) + 1        # 総フレーム数
interval_ms = 1000 / fps                    # 1コマあたりのミリ秒
t = np.linspace(start_time, end_time, steps)
t_pre = t[t < event_time]
t_post = t[t >= event_time]

stress = np.where(t - event_time >= 0, stress_i, 0)

# solution of ODE
s0 = stress_i   # ODEの引数として入れるためにこの形で定義
e0 = s0/E1      # ステップ応力を加えた直後の歪みはspring1の歪みのみであるため、初期条件として定義
sol = odeint(SLS1_stepStress, e0, t_post, args=(s0,infMod,tau))
strain_pre = np.zeros_like(t_pre)   # ステップ前の歪みはゼロ
strain_post = sol[:, 0]             # 歪み履歴
strain = np.concatenate([strain_pre, strain_post])

# 描画のためのスケーリング
e = strain/1.0          # 描画のためのスケーリング
s = stress/10**6        # 描画のためのスケーリング ([MPa]単位に変換)
e1 = s*10**6/E1         # [] バネ1の歪み
e2 = e - e1             # [] フォークト要素の歪み
s1 = E2*e2/10**6        # [MPa] フォークト要素のバネ2の応力
de2dt = np.gradient(e2, t)  # numpyを使ったフォークト要素の歪みの微分
s2 = eta*de2dt/10**6    # フォークト要素のダッシュポットの応力 ([MPa]単位に変換)

# グラフの初期設定
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

# --- 上段：応力 (Input) ---
ax1.set_xlim(start_time, end_time)
ax1.set_ylim(-np.max(s)*0.5, np.max(s)*1.8) # 縦軸を固定
ax1.set_ylabel('Applied stress, $\sigma$ /MPa')
ax1.set_title("SLS I model: step stress")
ax1.grid(True, ls='--')

line_stress, = ax1.plot([], [], color='red', lw=2, label='Step stress (t = {0:.1f} s)'.format(event_time))
line_stress_s1, = ax1.plot([], [], color='green', ls="dashed", lw=1, label='$\sigma$ (spring 2 in Voigt)')
line_stress_s2, = ax1.plot([], [], color='orange', ls="dashed", lw=1, label='$\sigma$ (dashpot in Voigt)')
ax1.legend(loc='upper right')

# --- 下段：歪み (Response) ---
ax2.set_xlim(start_time, end_time)
ax2.set_ylim(-np.max(e)*0.5, np.max(e)*1.5) # 縦軸を固定
ax2.set_xlabel('$t$ /s')
ax2.set_ylabel('Strain, $\epsilon$ /')
ax2.grid(True, ls='--')

line_strain, = ax2.plot([], [], color='blue', lw=2, label='Response to step strain')
line_strain_e1, = ax2.plot([], [], color='cyan', ls="dashed", lw=1, label='$\epsilon$ (spring 1)')
line_strain_e2, = ax2.plot([], [], color='gray', ls="dashed", lw=1, label='$\epsilon$ (Voigt)')
ax2.legend(loc='upper right')

# テキスト描画
var_text = r'$\sigma_0$ = {0:.2f} MPa, $E_1$ = {1:.1f} MPa, $E_2$ = {2:.1f} MPa, $\eta$ = {3:.1f} kPa s'.format(s0/10**6,E1/10**6,E2/10**6,eta/10**3)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'd$\epsilon$/d$t$ = ($\sigma_0$/$E_{{{\infty}}}$ - $\epsilon$)/$\tau$'
ax2.text(0.1, 0.9, eq_text, transform=ax2.transAxes)
res_text = r'$E_i$ = {0:.2f} MPa, $E_\infty$ = {1:.2f} MPa, $\tau$ = {2:.2f} s'.format(insMod/10**6, infMod/10**6, tau)
ax2.text(0.1, 0.8, res_text, transform=ax2.transAxes)

# アニメーション更新関数
def animate(i):
    # 応力データの更新
    line_stress.set_data(t[:i], s[:i])
    line_stress_s1.set_data(t[:i], s1[:i])
    line_stress_s2.set_data(t[:i], s2[:i])

    # 歪みデータの更新
    line_strain.set_data(t[:i], e[:i])
    line_strain_e1.set_data(t[:i], e1[:i])
    line_strain_e2.set_data(t[:i], e2[:i])

    return line_stress, line_stress_s1, line_stress_s2, line_strain, line_strain_e1, line_strain_e2

# アニメーション実行   
ani = animation.FuncAnimation(fig, animate, frames=steps, interval=interval_ms, blit=True)

savefile = './mp4/SLS1_stepStress.mp4'
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '30'])

plt.tight_layout()
plt.show()