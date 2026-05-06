# マクスウェルモデルの常微分方程式（振動歪み）
# タイムラインチャート版

import numpy as np
from scipy.integrate import odeint
from scipy import integrate
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib.animation as animation

'''
テキストの式(2.11)をベースに組み立てる
'''

# 変数の設定
try:
    E = float(input('modulus [MPa] (default = 0.4 MPa): '))*10**6
except ValueError:
    E = 4*10**5                 # [Pa] 弾性率
try:
    eta = float(input('viscosity [kPa s] (default = 100.0 kPa s): '))*10**3
except ValueError:
    eta = 10**5                 # [Pa s] 粘度

tau = eta/E

# 振動歪みの設定
try:
    log_freq_min = float(input('log(freq_min) for forced oscillation (default=-1.0): '))
except ValueError:
    log_freq_min = -1.0
try:
    log_freq_max = float(input('log(freq_max) for forced oscillation (default=0.5): '))
except ValueError:
    log_freq_max = 0.5
try:
    num_freq = int(input('number of frequency (default=5): '))
except ValueError:
    num_freq = 5
try:
    eamp = float(input('amplitude for sinusoidal strain [] (default=0.2): '))
except ValueError:
    eamp = 0.2

freq_list = np.logspace(log_freq_min, log_freq_max, num_freq)

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
t_start_pre = -2.0          # [s] 開始時間
t_event = 0.0               # [s] 振動歪みを加える時刻
t_duration_pre = t_event - t_start_pre          # [s] ステップ前の継続時間
fps = 60                    # 1秒あたりのフレーム数、30だと足りないので60に変更
interval_ms = 1000 / fps    # 1コマあたりのミリ秒
steps_pre = int(t_duration_pre * fps) + 1       # 総フレーム数
t_pre = np.linspace(t_start_pre, t_event, steps_pre)
strain_pre = np.zeros_like(t_pre)               # ステップ前の歪みはゼロ
stress_pre = np.zeros_like(t_pre)               # ステップ前の応力はゼロ
i_stress_pre = np.zeros_like(t_pre)
af_pre = np.zeros_like(t_pre) 
samp_pre = np.zeros_like(t_pre) 
pdiff_pre = np.zeros_like(t_pre)
t_start = t_pre[-1]

s0 = 0                  # ODEの初期条件として定義
is0 = 0                 # 応力の積分の初期条件
t_ani = t_pre
strain = strain_pre     # 周波数掃引の全ての入力信号（歪み）を格納
stress = stress_pre     # 周波数掃引の全ての出力信号（応力）を格納
i_stress = i_stress_pre
af_ani = af_pre         # 入力角周波数を格納（アニメーション用）
samp_ani = samp_pre     # 出力振幅の最大値を格納（アニメーション用
pdiff_ani = pdiff_pre   # 出力信号の位相を格納（アニメーション用）
samp_list = []          # 各周波数での出力振幅の最大値を格納
pdiff_list = []         # 各周波数での出力信号の位相を格納

# 各周波数での計算
for freq in freq_list:
    # データ準備
    t_duration = np.where(4.0/freq > 30.0, 2.0/freq, 4.0/freq)    # [s] 継続時間
#    t_duration = 4.0/freq    # [s] 継続時間
    steps = int(t_duration * fps) + 1   # 総フレーム数
    t_end = t_start + t_duration
    t = np.linspace(t_start, t_end, steps)
    t_ani = np.concatenate([t_ani, t]) 
#    print(t[0],t[-1],len(t))
    af = 2*np.pi*freq
    # 各周波数でのアニメーション期間中に同じ数値をずっと表示させるため
    af_ani_f = af*np.ones_like(t)       # 入力周波数を格納（アニメーション用）
    af_ani = np.concatenate([af_ani, af_ani_f]) 
    strain_f = eamp*np.sin(af*(t - t_start))        # 入力信号
    strain = np.concatenate([strain, strain_f])
    # ODEの解析
    sol = odeint(Maxwell_sinuStrain, s0, t - t_start, args=(eamp,af,E,tau)) # ODEの解
    stress_f = sol[:, 0]            # 応力履歴
    stress = np.concatenate([stress, stress_f]) 
    s0 = stress[-1]                 # 次のODE計算のために初期条件s0を更新
    i_stress_f = integrate.cumulative_trapezoid(stress_f, t, initial=0)  # scipyを使った応力の積分
    # 上で「initial=is0」としたいのだが、「0 or None to the initial argument is officially deprecated」とのことで以下の対応に変更
    i_stress_f += is0               # 前の最後のi_stressにスムーズに繋ぐために全体からis0を減算
    i_stress = np.concatenate([i_stress, i_stress_f])
    is0 = i_stress[-1]              # 次の応力の積分の計算のために確保
    # 位相差の計算
    strain_latter = strain_f[int(0.4*len(strain_f)):]     # 後半部分を抽出（前半は過渡応答を含むから）
    stress_latter = stress_f[int(0.4*len(stress_f)):]     # 後半部分を抽出（前半は過渡応答を含むから）
    samp = np.max(stress_latter)
    samp_ani_f = samp*np.ones_like(t)         # 入力周波数を格納（アニメーション用）
    samp_ani = np.concatenate([samp_ani, samp_ani_f]) 
    samp_list.append(samp)
    ind = getNearestIndex2value(stress_latter,0)          # 出力信号が0になるindexを抽出           
    pdiff = (180/np.pi)*np.arcsin(np.abs(strain_latter[ind])/eamp)
    pdiff_ani_f = pdiff*np.ones_like(t)         # 入力周波数を格納（アニメーション用）
    pdiff_ani = np.concatenate([pdiff_ani, pdiff_ani_f]) 
    pdiff_list.append(pdiff)
#    print(samp, pdiff)
    t_start = t[-1]

# 描画のためのスケーリング
e = strain/1.0     # 描画のためのスケーリング_
s = stress/10**6   # 描画のためのスケーリング ([MPa]単位に変換)
samp_ani = samp_ani/10**6
samp_max = np.max(samp_list)/10**6
e_s = s/(E/10**6)   # バネの歪み
e_d = i_stress/eta   # ダッシュポットの歪み（本来はこちらが正しい）
#e_d = e - e_s      # 簡単にはこちらでも良い    

# グラフの初期設定
fig = plt.figure(figsize=(8,5), tight_layout=True)
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
title_text = "Maxwell model: sinusoidal strain (frequecy sweep)"
ax1.set_title(title_text)
ax1.set_ylim(-3*eamp, 3*eamp)
ax2.set_ylim(-2*samp_max, 2*samp_max)
ax1.set_axisbelow(True)
ax2.set_axisbelow(True)
ax1.grid(True)
ax2.grid(False)
ax1.set_xlabel('$t$ /s')
ax1.set_ylabel('strain, $\epsilon$')
ax2.set_ylabel('stress, $\sigma$ /MPa')

# テキスト描画
var_text = r'$\epsilon_{{amp}}$ = {0:.3f}, $E$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s'.format(eamp,E/10**6,eta/10**3)
ax1.text(0.1, 0.25, var_text, transform=ax1.transAxes)
eq_text = r'd$\sigma$/d$t$ = -$\sigma$/$\tau$ + $E$d$\epsilon$/d$t$'
ax1.text(0.1, 0.15, eq_text, transform=ax1.transAxes)
res_text = r'$\tau$ = {0:.3f} s'.format(tau)
ax1.text(0.1, 0.05, res_text, transform=ax1.transAxes)

strain, = ax1.plot([], [], 'b', animated=True, label='strain, $\epsilon$ (input)')
strain_s, = ax1.plot([], [], 'g--', animated=True, label='strain, $\epsilon$ (spring)')
strain_d, = ax1.plot([], [], 'y--', animated=True, label='strain, $\epsilon$ (dashpot)')
stress, = ax2.plot([], [], 'r', animated=True, label='stress, $\sigma$ (output)')
# ここでは[],[]としているが、下で実際の値を入れている

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2)

samp_template = r'$\sigma_{{amp}}$ = %.3f MPa'
samp_text = ax1.text(0.35, 0.9, '', transform=ax1.transAxes)
pdiff_template = r'$\theta$ = %.1f$\degree$'
pdiff_text = ax1.text(0.35, 0.8, '', transform=ax1.transAxes)
af_template = r'$\omega$ = %.2f s$^{{-1}}$'
af_text = ax1.text(0.1, 0.8, '', transform=ax1.transAxes)
aft_template = r'$\omega\tau$ = %.2f'
aft_text = ax1.text(0.1, 0.7, '', transform=ax1.transAxes)

time_template = '$t$ = %.1f s'
time_text = ax1.text(0.1, 0.9, '', transform=ax1.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

# アニメーション更新関数
def animate(i):
    window_width = 10.0
    current_time = t_ani[i]

    # データの更新
    strain.set_data(t_ani[:i], e[:i])
    strain_s.set_data(t_ani[:i], e_s[:i])
    strain_d.set_data(t_ani[:i], e_d[:i])
    stress.set_data(t_ani[:i], s[:i])

    samp_text.set_text(samp_template % samp_ani[i])
    pdiff_text.set_text(pdiff_template % pdiff_ani[i])
    af_text.set_text(af_template % af_ani[i])
    aft_text.set_text(aft_template % (af_ani[i] * tau))

    # スクロール
    if current_time > window_width:
        ax1.set_xlim(current_time - window_width, current_time)
    else:
        ax1.set_xlim(t_start_pre, window_width)

    # 時刻表示
    time_text.set_text(time_template % (i/fps + t_start_pre))

    return strain, strain_s, strain_d, stress, samp_text, pdiff_text, af_text, aft_text, time_text

# アニメーション実行   
ani = animation.FuncAnimation(fig, animate, frames=len(t_ani), interval=interval_ms, blit=False, repeat=False)

'''
軸も更新（スクロール）したい場合は、blitはFalseに設定
Trueだと高速だが、背景や軸が更新されない
'''

savefile = './mp4/Maxwell_sinuStrain_timeLineChart.mp4'
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '60'])

plt.tight_layout()
plt.show()