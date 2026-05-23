# SLS Iの常微分方程式（振動歪み）
# 周波数応答（動的弾性率）

'''
How to use
% python3 SLS1_sinuStrain_dynamicModuli.py args[1]
args: -log
"-log"をつけると縦軸をログスケールに変換
何もついていないか、間違えたものがついている時はリニアスケールで表示
'''

import sys
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

if len(sys.argv) == 1:
    axisoption = ""
else:
    axisoption = sys.argv[1]
if axisoption == "-log":
    pass
else:
    axisoption = ""

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
    log_freq_min = float(input('log(freq_min) for forced oscillation (default=-1.2): '))
except ValueError:
    log_freq_min = -1.2
try:
    log_freq_max = float(input('log(freq_max) for forced oscillation (default=1.5): '))
except ValueError:
    log_freq_max = 1.5
try:
    num_freq = int(input('number of frequency (default=31): '))
except ValueError:
    num_freq = 31
try:
    eamp = float(input('amplitude for sinusoidal strain [] (default=0.2): '))
except ValueError:
    eamp = 0.2

freq_list = np.logspace(log_freq_min, log_freq_max, num_freq)
aft_array = 2 * np.pi * freq_list * tau

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
t_start_pre = -2.0          # [s] 開始時間
t_event = 0.0               # [s] 振動歪みを加える時刻
t_duration_pre = t_event - t_start_pre          # [s] ステップ前の継続時間
fps = 1000                    # 1秒あたりのフレーム数
# ここではフレーム数としては意味がないが、stepsを決めるために使っている。細かくデータを取るために1000を代入している
interval_ms = 1000 / fps    # 1コマあたりのミリ秒
steps_pre = int(t_duration_pre * fps) + 1       # 総フレーム数
t_pre = np.linspace(t_start_pre, t_event, steps_pre)
strain_pre = np.zeros_like(t_pre)               # ステップ前の歪みはゼロ
stress_pre = np.zeros_like(t_pre)               # ステップ前の応力はゼロ
dedt_pre = np.zeros_like(t_pre)
af_pre = np.zeros_like(t_pre) 
eamp_pre = np.zeros_like(t_pre) 
pdiff_pre = np.zeros_like(t_pre)
t_start = t_pre[-1]

s0 = 0                  # ODEの初期条件として定義
strain = strain_pre     # 周波数掃引の全ての入力信号（歪み）を格納
stress = stress_pre     # 周波数掃引の全ての出力信号（応力）を格納
dedt = dedt_pre
samp_list = []          # 各周波数での出力振幅の最大値を格納
pdiff_list = []         # 各周波数での出力信号の位相を格納

# 各周波数での計算
for freq in freq_list:
    # データ準備
#    t_duration = np.where(4.0/freq > 30.0, 2.0/freq, 4.0/freq)    # [s] 継続時間
    t_duration = 20.0/freq    # [s] 継続時間
    steps = int(t_duration * fps) + 1   # 総フレーム数
    t_end = t_start + t_duration
    t = np.linspace(t_start, t_end, steps)
#    print(t[0],t[-1],len(t))
    af = 2*np.pi*freq
    strain_f = eamp*np.sin(af*(t - t_start))        # 入力信号
    strain = np.concatenate([strain, strain_f])
    # ODEの解析
    sol = odeint(SLS1_sinuStrain, s0, t - t_start, args=(eamp,af,insMod,k,tau)) # ODEの解
    stress_f = sol[:, 0]            # 応力履歴
    stress = np.concatenate([stress, stress_f]) 
    s0 = stress[-1]                 # 次のODE計算のために初期条件e0を更新
    # 位相差の計算
    stress_latter = stress_f[int(0.4*len(stress_f)):]     # 後半部分を抽出（前半は過渡応答を含むから）
    strain_latter = strain_f[int(0.4*len(strain_f)):]     # 後半部分を抽出（前半は過渡応答を含むから）
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

# E', E"の計算
'''
E' = samp*cos(theta)/eamp
E" = samp*sin(theta)/eame
'''

cos_pdiff = np.cos(np.radians(pdiff_array))
sin_pdiff = np.sin(np.radians(pdiff_array))
strMod = samp_array * cos_pdiff / eamp
losMod = samp_array * sin_pdiff / eamp

if axisoption == "-log":
    strMod = strMod*10**6
    losMod = losMod*10**6
else:
    strMod = strMod
    losMod = losMod

fig = plt.figure(figsize=(8,5), tight_layout=True)
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
ax2.grid(ls='dotted')
title_text = "SLS I model: sinusoidal strain (frequecy sweep)"
ax1.set_title(title_text)
ax1.set_axisbelow(True)
ax1.set_xscale('log')
ax1.set_xlabel(r'$\omega\tau$')

# テキスト描画
var_text = r'$\epsilon_{{amp}}$ = {0:.2f}, $E_1$ = {1:.1f} MPa, $E_2$ = {2:.1f} MPa, $\eta$ = {3:.1f} kPa s'.format(eamp,E1/10**6,E2/10**6,eta/10**3)
eq_text = r'd$\sigma$/d$t$ = ($E_i$$\epsilon$ + $E_i$$\tau$ d$\epsilon$/d$t$ - $k$$\sigma$)/$\tau$'
res_text = r'$\tau$ = {0:.1f} s'.format(tau)

if axisoption == "-log":
    ax1.set_ylim(10**3, 10**7)
    ax1.set_ylabel(r'storage modulus, $E^{{\prime}}$ /Pa')
    ax2.set_ylim(10**3, 10**7)
    ax2.set_ylabel(r'loss modulus, $E^{{\prime\prime}}$ /Pa')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax1.text(0.1, 0.3, var_text, transform=ax1.transAxes)
    ax1.text(0.1, 0.2, eq_text, transform=ax1.transAxes)
    ax2.text(0.1, 0.1, res_text, transform=ax2.transAxes)
else:
    ax1.set_ylim(-0.05*np.max(strMod), 1.2*np.max(strMod))
    ax1.set_ylabel(r'storage modulus, $E^{{\prime}}$ /MPa')
    ax2.set_ylim(-0.05*np.max(strMod), 1.2*np.max(strMod))
    ax2.set_ylabel(r'loss modulus, $E^{{\prime\prime}}$ /MPa')
    ax1.text(0.15, 0.8, var_text, transform=ax1.transAxes)
    ax1.text(0.15, 0.7, eq_text, transform=ax1.transAxes)
    ax2.text(0.15, 0.6, res_text, transform=ax2.transAxes)

ax1.plot(aft_array,strMod, 'ro-', label=r'$E^{{\prime}}$')
ax2.plot(aft_array,losMod, 'bo-', label=r'$E^{{\prime\prime}}$')

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2)

if axisoption == "-log":
    savefile = "./png/SLS1_sinuStrain_dynamicModuli(log)_(tau={0:.1f}s).png".format(tau)
else:
    savefile = "./png/SLS1_sinuStrain_dynamicModuli(linear)_(tau={0:.1f}s).png".format(tau)

fig.savefig(savefile, dpi=300)

plt.show()