# フォークトモデルの常微分方程式（振動応力）
# 周波数応答（動的コンプライアンス）

'''
How to use
% python3 Maxwell_sinuStrain_dynamicModuli.py args[1]
args: -log
"-log"をつけると縦軸をログスケールに変換
何もついていないか、間違えたものがついている時はリニアスケールで表示
'''

import sys
import numpy as np
from scipy.integrate import odeint
from scipy import integrate
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
テキストの式(2.23)をベースに組み立てる
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
    num_freq = int(input('number of frequency (default=5): '))
except ValueError:
    num_freq = 31
try:
    samp = float(input('amplitude for sinusoidal stress [MPa] (default=0.2): '))*10**6
except ValueError:
    samp = 0.2*10**6

freq_list = np.logspace(log_freq_min, log_freq_max, num_freq)
aft_array = 2 * np.pi * freq_list * tau

# ODE解析で用いる関数の定義
def Voigt_sinuStress(e, t, samp, af, E, tau):
# e: 歪み, s: 応力, E: 弾性率, tau: 遅延時間
# ここではsampとafを指定し、この中でsの関数を作り振動応力を実現
    s = samp*np.sin(af*t)           # 振動応力の関数
    dedt = (s/E - e)/tau            # (2.23)
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
    sol = odeint(Voigt_sinuStress, e0, t - t_start, args=(samp,af,E,tau)) # ODEの解
    strain_f = sol[:, 0]            # 応力履歴
    strain = np.concatenate([strain, strain_f]) 
    e0 = strain[-1]                 # 次のODE計算のために初期条件e0を更新
    dedt_f = np.gradient(strain_f, t)
    dedt = np.concatenate([dedt, dedt_f]) 
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
title_text = "Voigt model: sinusoidal stress (frequecy sweep)"
ax1.set_title(title_text)
ax1.set_axisbelow(True)
ax1.set_xscale('log')
ax1.set_xlabel(r'$\omega\tau$')
if axisoption == "-log":
    ax1.set_ylim(10**(-8), 10**(-5))
    ax1.set_ylabel(r'storage compliance, $J^{{\prime}}$ /Pa$^{{-1}}$')
    ax2.set_ylim(10**(-8), 10**(-5))
    ax2.set_ylabel(r'loss compliance, $J^{{\prime\prime}}$ /$^{{-1}}$')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
else:
    ax1.set_ylim(-0.05*np.max(strComp), 1.2*np.max(strComp))
    ax1.set_ylabel(r'storage compliance, $J^{{\prime}}$ /MPa$^{{-1}}$')
    ax2.set_ylim(-0.05*np.max(strComp), 1.2*np.max(strComp))
    ax2.set_ylabel(r'loss compliance, $J^{{\prime\prime}}$ /MPa$^{{-1}}$')

var_text = r'$\sigma_{{amp}}$ = {0:.2f} MPa, $E$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s'.format(samp/10**6,E/10**6,eta/10**3)
ax1.text(0.15, 0.3, var_text, transform=ax1.transAxes)
eq_text = r'd$\epsilon$/d$t$ = ($\sigma$/$E$ - $\epsilon$)/$\tau$'
ax1.text(0.15, 0.2, eq_text, transform=ax1.transAxes)
res_text = r'$\tau$ = {0:.2f} s'.format(tau)
ax1.text(0.15, 0.1, res_text, transform=ax1.transAxes)

ax1.plot(aft_array,strComp, 'ro-', label=r'$J^{{\prime}}$')
ax2.plot(aft_array,losComp, 'bo-', label=r'$J^{{\prime\prime}}$')

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2)

if axisoption == "-log":
    savefile = "./png/Voigt_sinuStress_dynamicCompliance(log)_(tau={0:.1f}s).png".format(tau)
else:
    savefile = "./png/Voigt_sinuStress_dynamicCompliance(linear)_(tau={0:.1f}s).png".format(tau)

fig.savefig(savefile, dpi=300)

plt.show()