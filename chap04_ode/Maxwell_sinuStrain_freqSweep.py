# マクスウェルモデルの常微分方程式（振動歪み）
# 周波数応答（振幅と位相）

import numpy as np
from scipy.integrate import odeint
from scipy import integrate
import matplotlib.pyplot as plt

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
    log_freq_min = float(input('log(freq_min) for forced oscillation (default=-1.8): '))
except ValueError:
    log_freq_min = -1.8
try:
    log_freq_max = float(input('log(freq_max) for forced oscillation (default=1.2): '))
except ValueError:
    log_freq_max = 1.2
try:
    num_freq = int(input('number of frequency (default=30): '))
except ValueError:
    num_freq = 31
try:
    eamp = float(input('amplitude for sinusoidal strain [] (default=0.2): '))
except ValueError:
    eamp = 0.2

freq_list = np.logspace(log_freq_min, log_freq_max, num_freq)
aft_array = 2 * np.pi * freq_list * tau

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
fps = 1000                  # 1秒あたりのフレーム数
# ここではフレーム数としては意味がないが、stepsを決めるために使っている。細かくデータを取るために1000を代入している
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
strain = strain_pre     # 周波数掃引の全ての入力信号（歪み）を格納
stress = stress_pre     # 周波数掃引の全ての出力信号（応力）を格納
i_stress = i_stress_pre
samp_ani = samp_pre     # 出力振幅の最大値を格納（アニメーション用
pdiff_ani = pdiff_pre   # 出力信号の位相を格納（アニメーション用）
samp_list = []          # 各周波数での出力振幅の最大値を格納
pdiff_list = []         # 各周波数での出力信号の位相を格納

# 各周波数での計算
for freq in freq_list:
    # データ準備
#    t_duration = np.where(4.0/freq > 30.0, 2.0/freq, 4.0/freq)    # [s] 継続時間（アニメーション用）
    t_duration = 20.0/freq              # [s] 継続時間
    steps = int(t_duration * fps) + 1   # 総フレーム数
    t_end = t_start + t_duration
    t = np.linspace(t_start, t_end, steps)
#    print(t[0],t[-1],len(t))
    af = 2*np.pi*freq
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
    samp_list.append(samp)
    ind = getNearestIndex2value(stress_latter,0)          # 出力信号が0になるindexを抽出           
    pdiff = (180/np.pi)*np.arcsin(np.abs(strain_latter[ind])/eamp)
    pdiff_list.append(pdiff)
#    print(samp, pdiff)
    t_start = t[-1]

# 描画のためのスケーリング
samp_array = np.array(samp_list)/10**6
samp_max = np.max(samp_array)
pdiff_array = np.array(pdiff_list)

fig = plt.figure(figsize=(8,5), tight_layout=True)
ax1 = fig.add_subplot(111)
#ax1.grid()
title_text = "Maxwell model: sinusoidal strain (frequecy sweep)"
ax1.set_title(title_text)
ax1.set_axisbelow(True)
ax1.set_xlabel(r'$\omega\tau$')
ax1.set_ylim(-0.1*samp_max,1.2*samp_max)
ax1.set_ylabel(r'$\sigma_{{amp}}$ /MPa')
ax2 = ax1.twinx()
ax2.grid(ls='dotted')
ax2.set_ylim(-10,100)
ax2.set_ylabel(r'$\theta$ /$\degree$')

var_text = r'$\epsilon_{{amp}}$ = {0:.2f}, $E$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s'.format(eamp,E/10**6,eta/10**3)
ax1.text(0.05, 0.5, var_text, transform=ax1.transAxes)
eq_text = r'd$\sigma$/d$t$ = -$\sigma$/$\tau$ + $E$d$\epsilon$/d$t$'
ax1.text(0.05, 0.4, eq_text, transform=ax1.transAxes)
res_text = r'$\tau$ = {0:.2f} s'.format(tau)
ax1.text(0.05, 0.3, res_text, transform=ax1.transAxes)

ax1.plot(aft_array,samp_array, 'ro-', label=r'$\sigma_{{amp}}$')
ax1.set_xscale('log')
ax2.plot(aft_array,pdiff_array, 'bo-', label=r'$\theta$')

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2)

savefile = "./png/Maxwell_sinuStrain_freqSweep_(tau={0:.2f}s,).png".format(tau)
fig.savefig(savefile, dpi=300)

plt.show()