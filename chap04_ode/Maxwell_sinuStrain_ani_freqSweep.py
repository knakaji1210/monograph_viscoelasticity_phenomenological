# マクスウェルモデルの常微分方程式（振動歪み）（アニメーション付き）
# 周波数掃引版

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
int_stress_pre = np.zeros_like(t_pre)
af_pre = np.zeros_like(t_pre) 
samp_pre = np.zeros_like(t_pre) 
pdiff_pre = np.zeros_like(t_pre)
t_start = t_pre[-1]

s0 = 0                  # 初期条件として定義
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
e_d = i_stress/eta    # ダッシュポットの歪み
#e_d = e - e_s      # 簡単にはこちらでも良い
l = 0.1             # [m] 自然長
w = 0.5             # 要素の長さと幅の比率
# 直列なので自然長l0=2lとなっている
el = e*2*l          # [m] 全体の伸び
el_s = e_s*2*l      # [m] バネの伸び
el_d = e_d*2*l      # [m] ダッシュポットの伸び

# グラフの初期設定
fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111, xlabel='$t$ /s')
ax.grid()
title_text = "Maxwell model: sinusoidal strain (frequecy sweep)"
ax.set_title(title_text)
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
ax.set_xlim(-0.05,0.4)
ax.set_ylim(-5,5)

# 枠組み描画の準備
y_0 = [0, 0]
ax.plot([0, 0.08*l],y_0, c='b')
ax.plot(0,0,'ro', markersize='10')
ax.plot([0,2*l],[-3,-3], c='g')
ax.plot([0,0],[-2.8,-3.2], c='g')
ax.plot([2*l,2*l],[-2.8,-3.2], c='g')

# ダッシュポット描画の準備
x_d1 = [0.08*l, 0.92*l]
y_d1 = [w, w]
y_d2 = [-w, -w]
ax.plot(x_d1,y_d1, c='b')
ax.plot(x_d1,y_d2, c='b')
ax.plot([0.08*l,0.08*l],[w,-w], c='b')
rect = patches.Rectangle(xy=(0.08*l, -w), width=0.83*l, height=2*w, facecolor='y')
ax.add_patch(rect)

# テキスト描画
var_text = r'$\epsilon_{{amp}}$ = {0:.2f}, $E$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s'.format(eamp,E/10**6,eta/10**3)
ax.text(0.4, 0.9, var_text, transform=ax.transAxes)
eq_text = r'd$\sigma$/d$t$ = -$\sigma$/$\tau$ + $E$d$\epsilon$/d$t$'
ax.text(0.4, 0.8, eq_text, transform=ax.transAxes)
res_text = r'$\tau$ = {0:.2f} s'.format(tau)
ax.text(0.4, 0.7, res_text, transform=ax.transAxes)
ax.text(0.35, 0.15, '$l_0$', transform=ax.transAxes)
ax.text(0.75, 0.38, '$\epsilon$ (input)', transform=ax.transAxes)
ax.text(0.75, 0.28, '$\sigma$ (output)', transform=ax.transAxes)

samp_template = r'$\sigma_{{amp}}$ = %.3f MPa'
samp_text = ax.text(0.75, 0.52, '', transform=ax.transAxes)
pdiff_template = r'$\theta$ = %.1f$\degree$'
pdiff_text = ax.text(0.75, 0.45, '', transform=ax.transAxes)
af_template = r'$\omega$ = %.3f s$^{{-1}}$'
af_text = ax.text(0.1, 0.8, '', transform=ax.transAxes)
aft_template = r'$\omega\tau$ = %.2f'
aft_text = ax.text(0.1, 0.7, '', transform=ax.transAxes)

# 枠組みの描画
rod, = ax.plot([],[], 'b', animated=True)
point, = ax.plot([],[], 'ro', markersize='10', animated=True)

# バネの描画
rod_sp, = ax.plot([],[], 'b', animated=True)
triangle, = ax.plot([],[], 'b', animated=True)

# ダッシュポットの描画
rod_da, = ax.plot([],[], 'b', animated=True)
damper, = ax.plot([],[], 'b', lw=4, animated=True)
# ここでは[],[]としているが、下で***.set_dataで実際の値を入れている

# 歪みの描画
strain, = ax.plot([],[], 'b', lw=2, animated=True)
arrow_e_p, = ax.plot([],[], 'b', marker=9, markersize='10', animated=True)
arrow_e_n, = ax.plot([],[], 'b', marker=8, markersize='10', animated=True)

# 応力の描画
stress, = ax.plot([],[], 'r', lw=2, animated=True)
arrow_s_p, = ax.plot([],[], 'r', marker=9, markersize='10', animated=True)
arrow_s_n, = ax.plot([],[], 'r', marker=8, markersize='10', animated=True)

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return rod, point, triangle, rod_da, damper, strain, arrow_e_p, arrow_e_n, arrow_s_p, arrow_s_n, samp_text, pdiff_text, af_text, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    # 枠組み
    x_o = 7*l/4 + el[i]
    x_rod = [x_o, l/4+x_o]
    rod.set_data(x_rod,y_0)
    point.set_data([l/4+x_o],[0])
    # バネ
    x_sp = x_o-el_s[i]-l/2
    x_tri = np.linspace(x_sp, x_o,100)
    y_tri = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri - x_sp)/(x_o-x_sp)-np.pi/2+0.1))-1)
    triangle.set_data(x_tri,y_tri)
    # ダッシュポット
    x_da = l/2 + el_d[i]
    x_rod_da = [x_da, x_da+3*l/4]
    x_damp = x_da
    y_damp = 0.7*w
    x_damper = [x_damp, x_damp]
    y_damper = [y_damp, -y_damp]
    rod_da.set_data(x_rod_da,y_0)
    damper.set_data(x_damper,y_damper)
    # 歪み
    x_strain = [2*l, 2*l + el[i]]
    strain.set_data(x_strain,[-1,-1])
    if e[i] > 0:
        arrow_e_p.set_data([2*l + el[i]],[-1])
    else:
        arrow_e_n.set_data([2*l + el[i]],[-1])
    # 応力
    a = 1.0    # 見かけ上の振幅
    x_stress = [2*l, 2*l + a*s[i]]
    stress.set_data(x_stress,[-2,-2])
    if s[i] > 0:
        arrow_s_p.set_data([2*l + a*s[i]],[-2])
    else:
        arrow_s_n.set_data([2*l + a*s[i]],[-2])
    # パラメータ
    samp_text.set_text(samp_template % samp_ani[i])
    pdiff_text.set_text(pdiff_template % pdiff_ani[i])
    af_text.set_text(af_template % af_ani[i])
    aft_text.set_text(aft_template % (af_ani[i] * tau))
    # 時刻
    time_text.set_text(time_template % (i/fps + t_start_pre))
    return rod, point, triangle, rod_da, damper, strain, arrow_e_p, arrow_e_n, arrow_s_p, arrow_s_n, samp_text, pdiff_text, af_text, time_text

'''
y_triの中の重要部分は
x_tri1 = np.linspace(a, b,100)
のとき
(xtri - a)/(b - a)
になる 
'''

# アニメーション実行  
ani = animation.FuncAnimation(fig, update, frames=len(t_ani), 
                    init_func=init, blit=False, interval=interval_ms, repeat=False)

savefile = './mp4/Maxwell_sinuStrain_ani_freqSweep.mp4'
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '60'])

plt.tight_layout()
plt.show()
