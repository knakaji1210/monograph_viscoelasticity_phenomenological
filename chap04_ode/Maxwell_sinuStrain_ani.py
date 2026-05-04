# マクスウェルモデルの常微分方程式（振動歪み）（アニメーション付き）

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
e = strain/1.0      # 描画のためのスケーリング
s = stress/10**6    # 描画のためのスケーリング ([MPa]単位に変換)
s_max = stress_max/10**6
e_s = s/(E/10**6)   # バネの歪み
e_d = integral_stress/eta   # ダッシュポットの歪み
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
title_text = "Maxwell model: sinusoidal strain"
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
var_text = r'$\epsilon_{{amp}}$ = {0:.2f}, $f$ = {1:.3f} Hz, $E$ = {2:.1f} MPa, $\eta$ = {3:.1f} kPa s'.format(eamp,freq,E/10**6,eta/10**3)
ax.text(0.4, 0.9, var_text, transform=ax.transAxes)
eq_text = r'd$\sigma$/d$t$ = -$\sigma$/$\tau$ + $E$d$\epsilon$/d$t$'
ax.text(0.4, 0.8, eq_text, transform=ax.transAxes)
res_text = r'$\tau$ = {0:.1f} s, $\omega\tau$ = {1:.3f}'.format(tau, af*tau)
ax.text(0.4, 0.7, res_text, transform=ax.transAxes)
ax.text(0.35, 0.15, '$l_0$', transform=ax.transAxes)
ax.text(0.75, 0.38, '$\epsilon$ (input)', transform=ax.transAxes)
ax.text(0.75, 0.28, '$\sigma$ (output)', transform=ax.transAxes)
ax.text(0.1, 0.8, '$\omega$ = {0:.2f} s$^{{-1}}$'.format(af), transform=ax.transAxes)
samp_text = r'$\sigma_{{amp}}$ = {0:.3f} MPa'.format(s_max)
ax.text(0.75, 0.52, samp_text, transform=ax.transAxes)
phase_diff_text = r'$\theta$ = {0:.1f} $\degree$'.format(phase_diff)
ax.text(0.75, 0.45, phase_diff_text, transform=ax.transAxes)

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

# アニメーション更新関数
def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return rod, point, triangle, rod_da, damper, strain, arrow_e_p, arrow_e_n, time_text

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
    if el[i] > 0:
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
    # 時刻
    time_text.set_text(time_template % (t[i]))
    return rod, point, triangle, rod_da, damper, strain, arrow_e_p, arrow_e_n, time_text

'''
y_triの中の重要部分は
x_tri1 = np.linspace(a, b,100)
のとき
(xtri - a)/(b - a)
になる 
'''

# アニメーション実行  
ani = animation.FuncAnimation(fig, update, frames=steps, 
                    init_func=init, blit=True, interval=interval_ms, repeat=False)

savefile = './mp4/Maxwell_sinuStrain_ani_(f={0:.2f}Hz).mp4'.format(freq)
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '60'])

plt.tight_layout()
plt.show()