# SLS Iモデルの常微分方程式（振動歪み）

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation

'''
テキストの式(5.11)をベースに組み立てる
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
    eta = float(input('viscosity [kPa s] (default = 100.0 kPa s): '))*10**3
except ValueError:
    eta = 10**5               # [Pa s] 粘度

insMod = E1                 # [Pa] 瞬間弾性率
infMod = E1*E2/(E1+E2)      # [Pa] 緩和弾性率
k = insMod/infMod
tau = eta/E2                # [s] 緩和時間

# 振動歪みの設定
try:
    eamp = float(input('amplitude for sinusoidal strain [] (default=0.1): '))
except ValueError:
    eamp = 0.1
try:
    freq = float(input('frequency for sinusoidal strain [Hz] (default=1.91 Hz): '))
except ValueError:
    freq = 6/np.pi

af = 2*np.pi*freq

# ODE解析で用いる関数の定義
def SLS1_sinuStrain(s, t, eamp, af, insMod, k, tau):
# e: 歪み, s: 応力, insMod: 瞬間弾性率, tau: 緩和時間
# ここではeampとafを指定し、この中でeの関数を作り振動歪みを実現
    e = eamp*np.sin(af*t)
    dedt = eamp*af*np.cos(af*t)
    dsdt = (insMod*e + insMod*tau*dedt - k*s)/tau   # (5.11')
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
fps = 500            # 1秒あたりのフレーム数、30だと足りないので60に変更
# 高周波の場合はfps=500とかしないといけない
steps = int(time_duration * fps) + 1        # 総フレーム数
interval_ms = 1000 / fps                    # 1コマあたりのミリ秒
t = np.linspace(start_time, end_time, steps)
t_pre = t[t < event_time]
t_post = t[t >= event_time]

strain_post = eamp*np.sin(af*t)  # 振動歪みの関数
strain = np.where(t - event_time >= 0, strain_post, 0)

# ODEの解析
s0 = 0                              # 初期条件として定義
sol = odeint(SLS1_sinuStrain, s0, t_post, args=(eamp,af,insMod,k,tau))
stress_pre = np.zeros_like(t_pre)   # ステップ前の応力はゼロ
stress_post = sol[:, 0]             # 応力履歴
stress = np.concatenate([stress_pre, stress_post])

# 位相差の計算
strain_latter = strain[int(0.4*len(strain)):]     # 後半部分を抽出（前半は過渡応答を含むから）
stress_latter = stress[int(0.4*len(stress)):]     # 後半部分を抽出（前半は過渡応答を含むから）
stress_max = np.max(stress_latter)
ind = getNearestIndex2value(stress_latter,0)      # 出力信号が0になるindexを抽出           
phase_diff = (180/np.pi)*np.arcsin(np.abs(strain_latter[ind])/eamp)

# 描画のためのスケーリング
e = strain/1.0          # 描画のためのスケーリング
s = stress/10**6        # 描画のためのスケーリング ([MPa]単位に変換)
s_max = stress_max/10**6
e1 = s*10**6/E1         # [] バネ1の歪み
e2 = e - e1             # [] フォークト要素の歪み
s1 = E2*e2/10**6        # [MPa] フォークト要素のバネ2の応力
de2dt = np.gradient(e2, t)  # numpyを使ったフォークト要素の歪みの微分
s2 = eta*de2dt/10**6    # フォークト要素のダッシュポットの応力 ([MPa]単位に変換)

# グラフの初期設定
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

# --- 上段：歪み (Input) ---
ax1.set_xlim(start_time, end_time)
ax1.set_ylim(-np.max(e)*1.5, np.max(e)*1.5) # 縦軸を固定
ax1.set_ylabel('Applied strain, $\epsilon$ /')
ax1.set_title("SLS I model: sinusoidal strain")
ax1.grid(True, ls='--')

line_strain, = ax1.plot([], [], color='blue', lw=2, label='Response to step strain')
line_strain_e1, = ax1.plot([], [], color='cyan', ls="dashed", lw=1, label='$\epsilon$ (spring 1)')
line_strain_e2, = ax1.plot([], [], color='gray', ls="dashed", lw=1, label='$\epsilon$ (Voigt)')
ax1.legend(loc='upper right')

# --- 下段：応力 (Response) ---
ax2.set_xlim(start_time, end_time)
ax2.set_ylim(-np.max(s)*1.5, np.max(s)*1.5) # 縦軸を固定
ax2.set_xlabel('$t$ /s')
ax2.set_ylabel('Stress, $\sigma$ /MPa')
ax2.grid(True, ls='--')

line_stress, = ax2.plot([], [], color='red', lw=2, label='Step stress (t = {0:.1f} s)'.format(event_time))
line_stress_s1, = ax2.plot([], [], color='green', ls="dashed", lw=1, label='$\sigma$ (spring 2 in Voigt)')
line_stress_s2, = ax2.plot([], [], color='orange', ls="dashed", lw=1, label='$\sigma$ (dashpot in Voigt)')
ax2.legend(loc='upper right')

# テキスト描画
var_text = r'$\epsilon_{{amp}}$ = {0:.2f}, $E_1$ = {1:.1f} MPa, $E_2$ = {2:.1f} MPa, $\eta$ = {3:.1f} kPa s'.format(eamp,E1/10**6,E2/10**6,eta/10**3)
ax1.text(0.1, 0.9, var_text, transform=ax1.transAxes)
eq_text = r'd$\sigma$/d$t$ = ($E_i$$\epsilon$ + $E_i$$\tau$ d$\epsilon$/d$t$ - $k$$\sigma$)/$\tau$'
ax2.text(0.1, 0.9, eq_text, transform=ax2.transAxes)
res_text = r'$E_i$ = {0:.2f} MPa, $E_\infty$ = {1:.2f} MPa, $\tau$ = {2:.2f} s, $\omega\tau$ = {3:.3f}'.format(insMod/10**6, infMod/10**6, tau, af*tau)
ax2.text(0.1, 0.8, res_text, transform=ax2.transAxes)
samp_text = r'$\sigma_{{amp}}$ = {0:.3f} MPa'.format(s_max)
ax2.text(0.1, 0.35, samp_text, transform=ax2.transAxes)
phase_diff_text = r'$\theta$ = {0:.1f} $\degree$'.format(phase_diff)
ax2.text(0.1, 0.25, phase_diff_text, transform=ax2.transAxes)

# アニメーション更新関数
def animate(i):
    # 歪みデータの更新
    line_strain.set_data(t[:i], e[:i])
    line_strain_e1.set_data(t[:i], e1[:i])
    line_strain_e2.set_data(t[:i], e2[:i])

    # 応力データの更新
    line_stress.set_data(t[:i], s[:i])
    line_stress_s1.set_data(t[:i], s1[:i])
    line_stress_s2.set_data(t[:i], s2[:i])


    return line_strain, line_strain_e1, line_strain_e2, line_stress, line_stress_s1, line_stress_s2

# アニメーション実行   
ani = animation.FuncAnimation(fig, animate, frames=steps, interval=interval_ms, blit=True)

savefile = './mp4/SLS1_sinuStrain_(f={0:.2f}Hz).mp4'.format(freq)
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '30'])

plt.tight_layout()
plt.show()