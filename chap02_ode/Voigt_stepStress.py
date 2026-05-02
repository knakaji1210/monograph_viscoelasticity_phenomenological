# ordinary differential equation of Voigt model (step stress)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation

'''
テキストの式(2.24)をベースに組み立てる
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
    stress_i = float(input('step stress [MPa] (default = 0.05 MPa): '))*10**6
except ValueError:
    stress_i = 0.05*10**6         # [Pa] step stress

# ODE解析で用いる関数の定義
def Voigt_stepStress(e, t, s, E, tau):
# e: strain, s: stress, E: modulus, tau: retardation time
# ここでは下でargsとしてs0=stress_iを入れてステップ応力を実現
    dedt = (s/E - e)/tau      # (2.24)
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
sol = odeint(Voigt_stepStress, e0, t_post, args=(s0,E,tau))
strain_pre = np.zeros_like(t_pre)  # ステップ前の歪みはゼロ
strain_post = sol[:, 0]     # 歪み履歴
strain = np.concatenate([strain_pre, strain_post])

#dt = t[1] - t[0]   # 時間刻み
#dedt = np.array([0.0]+[(strain[k+1]-strain[k])/(t[k+1]-t[k]) for k in range(len(strain)-1)])   # 簡易的な歪みの微分
dedt = np.gradient(strain, t)                                                                   # numpyを使った歪みの微分

# scaling for figure
e = strain/1.0          # 描画のためのスケーリング
s = stress/10**6        # 描画のためのスケーリング ([MPa]単位に変換)
s_s = E*e/10**6         # バネの応力 ([MPa]単位に変換)
s_d = eta*dedt/10**6    # ダッシュポットの応力 ([MPa]単位に変換)

# 2. グラフの初期設定
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

# --- 上段：応力 (Input) ---
ax1.set_xlim(start_time, end_time)
ax1.set_ylim(-np.max(s)*0.5, np.max(s)*1.8) # 縦軸を固定
ax1.set_ylabel('Applied stress, $\sigma$ /MPa')
ax1.set_title("Voigt model: step stress")
ax1.grid(True, ls='--')

line_stress, = ax1.plot([], [], color='red', lw=2, label='Step stress (t = {0:.1f} s)'.format(event_time))
line_stress_s, = ax1.plot([], [], color='green', ls="dashed", lw=1, label='$\sigma$ (spring)')
line_stress_d, = ax1.plot([], [], color='orange', ls="dashed", lw=1, label='$\sigma$ (dashpot)')
ax1.legend(loc='upper right')

# --- 下段：歪み (Response) ---
ax2.set_xlim(start_time, end_time)
ax2.set_ylim(-np.max(e)*0.5, np.max(e)*1.5) # 縦軸を固定
ax2.set_xlabel('$t$ /s')
ax2.set_ylabel('Strain, $\epsilon$ /')
ax2.grid(True, ls='--')

line_strain, = ax2.plot([], [], color='blue', lw=2, label='Response to step strain')
ax1.legend(loc='upper right')

var_text = r'$\sigma_0$ = {0:.2f} MPa, $E$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s'.format(stress_i/10**6,E/10**6,eta/10**3)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'd$\epsilon$/d$t$ = ($\sigma_0$/$E$ - $\epsilon$)/$\tau$'
ax2.text(0.1, 0.9, eq_text, transform=ax2.transAxes)
res_text = r'$\tau$ = {0:.1f} s'.format(tau)
ax2.text(0.1, 0.8, res_text, transform=ax2.transAxes)

# 3. アニメーション更新関数
def animate(i):
    # 応力データの更新
    line_stress.set_data(t[:i], s[:i])
    line_stress_s.set_data(t[:i], s_s[:i])
    line_stress_d.set_data(t[:i], s_d[:i])

    # 歪みデータの更新
    line_strain.set_data(t[:i], e[:i])

    return line_stress, line_stress_s, line_stress_d, line_strain

# アニメーション実行   
ani = animation.FuncAnimation(fig, animate, frames=steps, interval=interval_ms, blit=True)

savefile = './mp4/Voigt_stepStress.mp4'
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '30'])

plt.tight_layout()
plt.show()