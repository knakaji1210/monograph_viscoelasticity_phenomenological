# SLS IIモデルの常微分方程式（振動応力）
# 周波数応答（動的コンプライアンス）

'''
How to use
% python3 SLS1_sinuStress_dynamicComp.py args[1]
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
テキストの式(5.39)をベースに組み立てる
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
    eta = 10**5             # [Pa s] 粘度

insMod = E1+E2              # [Pa] 瞬間弾性率
infMod = E2                 # [Pa] 緩和弾性率
k = insMod/infMod
tau = eta/E1                # [s] 緩和時間

# 振動応力の設定
try:
    log_freq_min = float(input('log(freq_min) for forced oscillation (default=-1.8): '))
except ValueError:
    log_freq_min = -1.8
try:
    log_freq_max = float(input('log(freq_max) for forced oscillation (default=1.2): '))
except ValueError:
    log_freq_max = 1.2
try:
    num_freq = int(input('number of frequency (default=31): '))
except ValueError:
    num_freq = 31
try:
    samp = float(input('amplitude for sinusoidal stress [MPa] (default=0.2): '))*10**6
except ValueError:
    samp = 0.2*10**6

freq_list = np.logspace(log_freq_min, log_freq_max, num_freq)
aft_array = 2 * np.pi * freq_list * tau

# ODE解析で用いる関数の定義
def SLS2_sinuStress(e, t, samp, af, insMod, k, tau):
# e: 歪み, s: 応力, insMod: 瞬間弾性率, tau: 緩和時間
# ここではsampとafを指定し、この中でsの関数を作り振動応力を実現
    s = samp * np.sin(af * t)
    dsdt = samp*af*np.cos(af*t)
    dedt = (s/insMod +tau*dsdt/insMod - e/k)/tau  # (5.39)
    return dedt

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

e0 = 0                  # ODEの初期条件として定義
strain = strain_pre     # 周波数掃引の全ての入力信号（歪み）を格納
stress = stress_pre     # 周波数掃引の全ての出力信号（応力）を格納
dedt = dedt_pre
eamp_list = []          # 各周波数での出力振幅の最大値を格納
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
    stress_f = samp*np.sin(af*(t - t_start))        # 入力信号
    stress = np.concatenate([stress, stress_f])
    # ODEの解析
    sol = odeint(SLS2_sinuStress, e0, t - t_start, args=(samp,af,insMod,k,tau)) # ODEの解
    strain_f = sol[:, 0]            # 応力履歴
    strain = np.concatenate([strain, strain_f]) 
    e0 = strain[-1]                 # 次のODE計算のために初期条件e0を更新
    # 位相差の計算
    strain_latter = strain_f[int(0.4*len(strain_f)):]     # 後半部分を抽出（前半は過渡応答を含むから）
    stress_latter = stress_f[int(0.4*len(stress_f)):]     # 後半部分を抽出（前半は過渡応答を含むから）
    eamp = np.max(strain_latter)
    eamp_list.append(eamp)
    ind = getNearestIndex2value(strain_latter,0)          # 出力信号が0になるindexを抽出           
    pdiff = (180/np.pi)*np.arcsin(np.abs(stress_latter[ind])/samp)
    pdiff_list.append(pdiff)
#    print(samp, pdiff)
    t_start = t[-1]

# 描画のためのスケーリング
eamp_array = np.array(eamp_list)/1.0
eamp_max = np.max(eamp_array)
pdiff_array = np.array(pdiff_list)

# J', J"の計算
'''
J' = eamp*cos(theta)/samp
J" = eamp*sin(theta)/samp
'''

cos_pdiff = np.cos(np.radians(pdiff_array))
sin_pdiff = np.sin(np.radians(pdiff_array))
strComp = eamp_array * cos_pdiff / samp
losComp = eamp_array * sin_pdiff / samp

if axisoption == "-log":
    strComp = strComp
    losComp = losComp
else:
    strComp = strComp * 10**6
    losComp = losComp * 10**6

fig = plt.figure(figsize=(8,5), tight_layout=True)
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
ax2.grid(ls='dotted')
title_text = "SLS II model: sinusoidal stress (frequecy sweep)"
ax1.set_title(title_text)
ax1.set_axisbelow(True)
ax1.set_xscale('log')
ax1.set_xlabel(r'$\omega\tau$')

# テキスト描画
var_text = r'$\sigma_{{amp}}$ = {0:.2f} MPa, $E_1$ = {1:.1f} MPa, $E_2$ = {2:.1f} MPa, $\eta$ = {3:.1f} kPa s'.format(samp/10**6,E1/10**6,E2/10**6,eta/10**3)
eq_text = r'd$\epsilon$/d$t$ = ($\sigma$/$E_i$ + $\tau$/$E_i$ d$\sigma$/d$t$ -$\epsilon$/$k$)/$\tau$'
res_text = r'$\tau$ = {0:.2f} s'.format(tau)

if axisoption == "-log":
    ax1.set_ylim(10**(-8), 10**(-5))
    ax1.set_ylabel(r'storage compliance, $J^{{\prime}}$ /Pa$^{{-1}}$')
    ax2.set_ylim(10**(-8), 10**(-5))
    ax2.set_ylabel(r'loss compliance, $J^{{\prime\prime}}$ /$^{{-1}}$')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax1.text(0.15, 0.3, var_text, transform=ax1.transAxes)
    ax1.text(0.15, 0.2, eq_text, transform=ax1.transAxes)
    ax2.text(0.15, 0.1, res_text, transform=ax2.transAxes)
else:
    ax1.set_ylim(-0.05*np.max(strComp), 1.2*np.max(strComp))
    ax1.set_ylabel(r'storage compliance, $J^{{\prime}}$ /MPa$^{{-1}}$')
    ax2.set_ylim(-0.05*np.max(strComp), 1.2*np.max(strComp))
    ax2.set_ylabel(r'loss compliance, $J^{{\prime\prime}}$ /MPa$^{{-1}}$')
    ax1.text(0.4, 0.8, var_text, transform=ax1.transAxes)
    ax1.text(0.4, 0.7, eq_text, transform=ax1.transAxes)
    ax2.text(0.4, 0.6, res_text, transform=ax2.transAxes)

ax1.plot(aft_array,strComp, 'ro-', label=r'$J^{{\prime}}$')
ax2.plot(aft_array,losComp, 'bo-', label=r'$J^{{\prime\prime}}$')

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2)

if axisoption == "-log":
    savefile = "./png/SLS2_sinuStress_dynamicCompliance(log)_(tau={0:.1f}s).png".format(tau)
else:
    savefile = "./png/SLS2_sinuStress_dynamicCompliance(linear)_(tau={0:.1f}s).png".format(tau)

fig.savefig(savefile, dpi=300)

plt.show()