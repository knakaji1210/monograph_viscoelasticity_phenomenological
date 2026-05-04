# ダッシュポット要素の常微分方程式（振動応力）

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation

'''
テキストの式(2.5)をベースに組み立てる
'''

# 変数の設定
try:
    eta = float(input('viscosity [kPa s] (default = 100.0 kPa s): '))*10**3
except ValueError:
    eta = 10**5             # [Pa s] 粘度

# 振動応力の設定
try:
    samp = float(input('amplitude for sinusoidal stress [MPa] (default=0.05): '))*10**6
except ValueError:
    samp = 0.05*10**6
try:
    freq = float(input('frequency for sinusoidal stress [Hz] (default=0.5 Hz): '))
except ValueError:
    freq = 0.5

af = 2*np.pi*freq

# ODE解析で用いる関数の定義
def dashpot_sinuStress(e, t, samp, af, eta):
# e: 歪み, s: 応力, eta: 粘度
# ここではsampとafを指定し、この中でsの関数を作り振動応力を実現
    s = samp*np.sin(af*t)           # 振動応力の関数
    dedt = s/eta                    # (1.5)
    return dedt

# データ準備
start_time = -2.0   # 開始時間
end_time = 4/freq   # 終了時間
event_time = 0.0    # 振動応力を加える時刻
time_duration = end_time - start_time       # [s] 継続時間
time_duration_pre = event_time - start_time # [s] 振動前の継続時間
time_duration_post = end_time - event_time  # [s] 振動後の継続時間
fps = 60            # 1秒あたりのフレーム数、30だと足りないので60に変更
steps = int(time_duration * fps) + 1        # 総フレーム数
interval_ms = 1000 / fps                    # 1コマあたりのミリ秒
t = np.linspace(start_time, end_time, steps)
t_pre = t[t < event_time]
t_post = t[t >= event_time]

stress_post = samp*np.sin(af*t)  # 振動応力の関数
stress = np.where(t - event_time >= 0, stress_post, 0)

# ODEの解析
e0 = 0.0                            # 初期条件として定義
sol = odeint(dashpot_sinuStress, e0, t_post, args=(samp,af,eta)) # ODEの解
strain_pre = np.zeros_like(t_pre)   # ステップ前の歪みはゼロ
strain_post = sol[:, 0]             # 歪み履歴
strain = np.concatenate([strain_pre, strain_post])

# 描画のためのスケーリング
e = strain/1.0     # 描画のためのスケーリング
s = stress/10**6   # 描画のためのスケーリング ([MPa]単位に変換)

# グラフの初期設定
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

# --- 上段：応力 (Input) ---
ax1.set_xlim(start_time, end_time)
ax1.set_ylim(-np.max(s)*1.5, np.max(s)*1.5) # 縦軸を固定
ax1.set_ylabel('Applied stress, $\sigma$ /MPa')
ax1.set_title("dashpot (Newton's viscosity): sinusoidal stress")
ax1.grid(True, ls='--')

line_stress, = ax1.plot([], [], color='red', lw=2, label='Sinusoidal stress (t = {0:.1f} s)'.format(event_time))
ax1.legend(loc='upper right')

# --- 下段：歪み (Response) ---
ax2.set_xlim(start_time, end_time)
ax2.set_ylim(-np.max(e)*1.5, np.max(e)*1.5) # 縦軸を固定
ax2.set_xlabel('$t$ /s')
ax2.set_ylabel('Strain, $\epsilon$ /')
ax2.grid(True, ls='--')

line_strain, = ax2.plot([], [], color='blue', lw=2, label='Response to sinusoidal stress')
ax2.legend(loc='upper right')

# テキスト描画
var_text = r'$\sigma_{{amp}}$ = {0:.2f} MPa, $f$ = {1:.1f} Hz, $\eta$ = {2:.1f} kPa s'.format(samp/10**6,freq,eta/10**3)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'd$\epsilon$/d$t$ = $\sigma$/$\eta$'
ax2.text(0.1, 0.9, eq_text, transform=ax2.transAxes)

# アニメーション更新関数
def animate(i):
    # 応力データの更新
    line_stress.set_data(t[:i], s[:i])
    
    # 歪みデータの更新
    line_strain.set_data(t[:i], e[:i])
    return line_stress, line_strain

# アニメーション実行   
ani = animation.FuncAnimation(fig, animate, frames=steps, interval=interval_ms, blit=True)

savefile = './mp4/dashpot_sinuStress_(f={0:.2f}Hz).mp4'.format(freq)
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '60'])

plt.tight_layout()
plt.show()