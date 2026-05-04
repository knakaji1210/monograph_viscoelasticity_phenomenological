# ダッシュポット要素の常微分方程式（振動応力）（アニメーション付き）

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib.animation as animation

'''
テキストの式(1.5)をベースに組み立てる
'''

# 変数の設定
try:
    eta = float(input('viscosity [kPa s] (default = 100.0 kPa s): '))*10**3
except ValueError:
    eta = 10**5           # [Pa s] 粘度

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

def getNearestIndex2value(list,value):
    index = np.abs(np.array(list) -value).argsort()[0].tolist()
    return index

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

# 位相差の計算
strain_latter = strain[int(0.4*len(strain)):]     # 後半部分を抽出（前半は過渡応答を含むから）
stress_latter = stress[int(0.4*len(stress)):]     # 後半部分を抽出（前半は過渡応答を含むから）
strain_max = np.max(strain_latter)
ind = getNearestIndex2value(strain_latter,strain_max/2)      # 出力信号が振幅の1/2になるindexを抽出           
phase_diff = (180/np.pi)*np.arcsin(np.abs(stress_latter[ind])/samp)

# 描画のためのスケーリング
e = strain/1.0     # 描画のためのスケーリング
e_max = strain_max/1.0
s = stress/10**6   # 描画のためのスケーリング ([MPa]単位に変換)
l = 0.1            # [m] 自然長
w = 0.5            # 要素の長さと幅の比率
el = e*l           # [m] 全体の伸び

# グラフの初期設定
fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111)
ax.grid()
title_text = "dashpot (Newton's viscosity): sinusoidal stress"
ax.set_title(title_text)
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
ax.set_xlim(-0.05,0.3)
ax.set_ylim(-5,5)

# 枠組みの描画
y_0 = [0, 0]
ax.plot([0, 0.08*l],y_0, c='b')
ax.plot(0,0, 'ro', markersize='10')

# ダッシュポット描画の準備
x_d1 = [0.08*l, 0.92*l]
y_d1 = [w, w]
y_d2 = [-w, -w]
ax.plot(x_d1,y_d1, c='b')
ax.plot(x_d1,y_d2, c='b')
ax.plot([0.08*l,0.08*l],[w,-w], c='b')
rect = patches.Rectangle(xy=(0.08*l, -w), width=0.83*l, height=2*w, facecolor='y')
ax.add_patch(rect)
ax.plot([0,l],[-3,-3], c='g')
ax.plot([0,0],[-2.8,-3.2], c='g')
ax.plot([l,l],[-2.8,-3.2], c='g')

# テキスト描画
var_text = r'$\sigma_{{amp}}$ = {0:.2f} MPa, $f$ = {1:.3f} Hz, $\eta$ = {2:.1f} kPa s'.format(samp/10**6,freq,eta/10**3)
ax.text(0.4, 0.9, var_text, transform=ax.transAxes)
eq_text = r'd$\epsilon$/d$t$ = $\sigma$/$\eta$'
ax.text(0.4, 0.8, eq_text, transform=ax.transAxes)
ax.text(0.3, 0.15, '$l_0$', transform=ax.transAxes)
ax.text(0.6, 0.28, '$\epsilon$ (output)', transform=ax.transAxes)
ax.text(0.6, 0.38, '$\sigma$ (input)', transform=ax.transAxes)
eamp_text = r'$\epsilon_{{amp}}$ = {0:.3f}'.format(e_max)
ax.text(0.75, 0.52, eamp_text, transform=ax.transAxes)
phase_diff_text = r'$\theta$ = {0:.1f} $\degree$'.format(phase_diff)
ax.text(0.75, 0.45, phase_diff_text, transform=ax.transAxes)

# ダッシュポットの描画
rod, = ax.plot([],[], 'b', animated=True)
damper, = ax.plot([],[], 'b', lw=4, animated=True)
point, = ax.plot([], [], 'ro', markersize='10', animated=True)
# ここでは[],[]としているが、下で***.set_dataで実際の値を入れている

# 応力の描画
stress_bar, = ax.plot([],[], 'r', lw=2, animated=True)
arrow_s_p, = ax.plot([],[], 'r', marker=9, markersize='10', animated=True)
arrow_s_n, = ax.plot([],[], 'r', marker=8, markersize='10', animated=True)

# 歪みの描画
strain_bar, = ax.plot([],[], 'b', lw=2, animated=True)
arrow_e_p, = ax.plot([],[], 'b', marker=9, markersize='10', animated=True)
arrow_e_n, = ax.plot([],[], 'b', marker=8, markersize='10', animated=True)

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

# アニメーション更新関数
def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return rod, damper, point, stress_bar, arrow_s_p, arrow_s_n, strain_bar, arrow_e_p, arrow_e_n, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    # 枠組み
    x_rod = [l/2 + el[i], l + el[i]]
    # ダッシュポット
    x_damp = (l/2)+el[i]
    y_damp = 0.7*w
    x_damper = [x_damp, x_damp]
    y_damper = [y_damp, -y_damp]
    rod.set_data(x_rod,y_0)
    damper.set_data(x_damper,y_damper)
    point.set_data([l + el[i]],[0])
    # 応力
    a = 0.5    # 見かけ上の振幅補正
    x_stress = [l, l + a*s[i]]
    stress_bar.set_data(x_stress,[-1,-1])
    if s[i] > 0:
        arrow_s_p.set_data([l + a*s[i]],[-1])
    else:
        arrow_s_n.set_data([l + a*s[i]],[-1])
    # 歪み
    x_strain = [l, l + el[i]]
    strain_bar.set_data(x_strain,[-2,-2])
    if el[i] > 0:
        arrow_e_p.set_data([l + el[i]],[-2])
    else:
        arrow_e_n.set_data([l + el[i]],[-2])
    time_text.set_text(time_template % (t[i]))
    return rod, damper, point, stress_bar, arrow_s_p, arrow_s_n, strain_bar, arrow_e_p, arrow_e_n, time_text

# アニメーション実行  
ani = animation.FuncAnimation(fig, update, frames=steps, 
                    init_func=init, blit=True, interval=interval_ms, repeat=False)

savefile = './mp4/dashpot_sinuStress_ani_(f={0:.2f}Hz).mp4'.format(freq)
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '60'])

plt.tight_layout()
plt.show()