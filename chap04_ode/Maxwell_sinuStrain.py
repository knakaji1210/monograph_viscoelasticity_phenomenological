# マクスウェルモデルの常微分方程式（振動歪み）

import numpy as np
from scipy.integrate import odeint
from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib.animation as animation

'''
テキストの式(2.11)をベースに組み立てる
'''

# 変数の設定
try:
    E = float(input('modulus [MPa] (default = 0.2 MPa): '))*10**6
except ValueError:
    E = 2*10**5                 # [Pa] 弾性率
try:
    eta = float(input('viscosity [kPa s] (default = 100.0 kPa s): '))*10**3
except ValueError:
    eta = 10**5                 # [Pa s] 粘度

tau = eta/E                     # [s] 緩和時間

# 振動歪みの設定
try:
    eamp = float(input('amplitude for sinusoidal strain [] (default=0.2): '))
except ValueError:
    eamp = 0.2
try:
    freq = float(input('frequency for sinusoidal strain [Hz] (default=0.3813 Hz): '))
except ValueError:
    freq = 1/np.pi

af = 2*np.pi*freq

# ODE解析で用いる関数の定義
def Maxwell_sinuStrain(s, t, eamp, af, E, tau):
# e: 歪み, s: 応力, E: 弾性率, tau: 緩和時間
# ここではeampとafを指定し、この中でeの関数を作り振動歪みを実現
    e = eamp*np.sin(af*t)           # 振動歪みの関数
    dedt = eamp*af*np.cos(af*t)     # 振動歪みの関数の時間微分
    dsdt = E*dedt-s/tau             # (1.11)
    return dsdt

def getNearestIndex2value(list,value):
    index = np.abs(np.array(list) -value).argsort()[0].tolist()
    return index

# データ準備
start_time = np.where(1.0/freq > 2.0, -2.0, -1.0/freq)      # 開始時間
#print('start_time = {0:.2f} s'.format(start_time))  
end_time = np.where(4.0/freq > 20.0, 2.0/freq, 4.0/freq)    # 終了時間
#print('end_time = {0:.2f} s'.format(end_time))
event_time = 0.0    # ステップ歪みを加える時刻
time_duration = end_time - start_time       # [s] 継続時間
time_duration_pre = event_time - start_time # [s] ステップ前の継続時間
time_duration_post = end_time - event_time  # [s] ステップ後の継続時間
fps = 60            # 1秒あたりのフレーム数、30だと足りないので60に変更
steps = int(time_duration * fps) + 1        # 総フレーム数
interval_ms = 1000 / fps                    # 1コマあたりのミリ秒
t = np.linspace(start_time, end_time, steps)
t_pre = t[t < event_time]
t_post = t[t >= event_time]

strain_post = eamp*np.sin(af*t)  # 振動歪みの関数
strain = np.where(t - event_time >= 0, strain_post, 0)

# ODEの解析
s0 = 0                              # 初期条件として定義
sol = odeint(Maxwell_sinuStrain, s0, t_post, args=(eamp,af,E,tau)) # ODEの解
stress_pre = np.zeros_like(t_pre)   # ステップ前の応力はゼロ
stress_post = sol[:, 0]             # 応力履歴
stress = np.concatenate([stress_pre, stress_post])

integral_stress = integrate.cumulative_trapezoid(stress, t, initial=0)  # scipyを使った応力の積分

# 位相差の計算
strain_latter = strain[int(0.4*len(strain)):]     # 後半部分を抽出（前半は過渡応答を含むから）
stress_latter = stress[int(0.4*len(stress)):]     # 後半部分を抽出（前半は過渡応答を含むから）
stress_max = np.max(stress_latter)
ind = getNearestIndex2value(stress_latter,0)      # 出力信号が0になるindexを抽出           
phase_diff = (180/np.pi)*np.arcsin(np.abs(strain_latter[ind])/eamp)

# 描画のためのスケーリング
e = strain/1.0     # 描画のためのスケーリング_
s = stress/10**6   # 描画のためのスケーリング ([MPa]単位に変換)
s_max = stress_max/10**6
e_s = s/(E/10**6)  # バネの歪み
e_d = integral_stress/eta   # ダッシュポットの歪み

# グラフの初期設定
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

# --- 上段：歪み (Input) ---
ax1.set_xlim(start_time, end_time)
ax1.set_ylim(-np.max(e)*1.5, np.max(e)*1.5) # 縦軸を固定
ax1.set_ylabel('Applied Strain, $\epsilon$ /')
ax1.set_title("Maxwell model: sinusoidal strain")
ax1.grid(True, ls='--')

line_strain, = ax1.plot([], [], color='blue', lw=2, label='Step strain (t = {0:.1f} s)'.format(event_time))
line_strain_s, = ax1.plot([], [], color='green', ls="dashed", lw=1, label='$\epsilon$ (spring)')
line_strain_d, = ax1.plot([], [], color='orange', ls="dashed", lw=1, label='$\epsilon$ (dashpot)')
ax1.legend(loc='upper right')

# --- 下段：応力 (Response) ---
ax2.set_xlim(start_time, end_time)
ax2.set_ylim(-np.max(s)*1.5, np.max(s)*1.5) # 縦軸を固定
ax2.set_xlabel('$t$ /s')
ax2.set_ylabel('Stress, $\sigma$ /MPa')
ax2.grid(True, ls='--')

line_stress, = ax2.plot([], [], color='red', lw=2, label='Response to step strain')
ax2.legend(loc='upper right')

# テキスト描画
var_text = r'$\epsilon_{{amp}}$ = {0:.2f}, $f$ = {1:.3f} Hz, $E$ = {2:.1f} MPa, $\eta$ = {3:.1f} kPa s'.format(eamp,freq,E/10**6,eta/10**3)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'd$\sigma$/d$t$ = -$\sigma$/$\tau$ + $E$d$\epsilon$/d$t$'
ax2.text(0.1, 0.9, eq_text, transform=ax2.transAxes)
res_text = r'$\tau$ = {0:.3f} s, $\omega\tau$ = {1:.3f}'.format(tau, af*tau)
ax2.text(0.1, 0.8, res_text, transform=ax2.transAxes)
samp_text = r'$\sigma_{{amp}}$ = {0:.3f} MPa'.format(s_max)
ax2.text(0.1, 0.35, samp_text, transform=ax2.transAxes)
phase_diff_text = r'$\theta$ = {0:.1f} $\degree$'.format(phase_diff)
ax2.text(0.1, 0.25, phase_diff_text, transform=ax2.transAxes)

# アニメーション更新関数
def animate(i):
    # 応力データの更新
    line_strain.set_data(t[:i], e[:i])
    line_strain_s.set_data(t[:i], e_s[:i])
    line_strain_d.set_data(t[:i], e_d[:i])

    # 歪みデータの更新
    line_stress.set_data(t[:i], s[:i])
    return line_strain, line_strain_s, line_strain_d, line_stress

# アニメーション実行   
ani = animation.FuncAnimation(fig, animate, frames=steps, interval=interval_ms, blit=True)

savefile = './mp4/Maxwell_sinuStrain_(f={0:.2f}Hz).mp4'.format(freq)
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '30'])

plt.tight_layout()
plt.show()