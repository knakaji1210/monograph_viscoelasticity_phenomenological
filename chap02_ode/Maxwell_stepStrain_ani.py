# マクスウェルモデルの常微分方程式（ステップ歪み）（アニメーション付き）

import numpy as np
from scipy.integrate import odeint
from scipy import integrate
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib.animation as animation

'''
テキストの式(2.12)をベースに組み立てる
'''

# 変数の設定
try:
    E = float(input('modulus [MPa] (default = 0.2 MPa): '))*10**6
except ValueError:
    E = 2*10**5                 # [Pa] 弾性率
try:
    eta = float(input('viscosity [kPa s] (default = 500.0 kPa s): '))*10**3
except ValueError:
    eta = 5*10**5               # [Pa s] 粘度

tau = eta/E                     # [s] 緩和時間

# 初期条件の設定
try:
    strain_i = float(input('step strain [] (default = 0.2): '))
except ValueError:
    strain_i = 0.2               # [] ステップ歪み

# ODE解析で用いる関数の定義
def Maxwell_stepStrain(s, t, tau):
# e: 歪み, s: 応力, tau: 緩和時間
# ステップ歪みe0=strain_iを加えるので、式(2.11)のde/dtの項が0となっている
    dsdt = -s/tau               # (2.12)
    return dsdt

# データ準備
start_time = -2.0   # 開始時間
end_time = 8.0      # 終了時間
event_time = 0.0    # ステップ歪みを加える時刻
time_duration = end_time - start_time       # [s] 継続時間
time_duration_pre = event_time - start_time # [s] ステップ前の継続時間
time_duration_post = end_time - event_time  # [s] ステップ後の継続時間
fps = 30            # 1秒あたりのフレーム数
steps = int(time_duration * fps) + 1        # 総フレーム数
interval_ms = 1000 / fps                    # 1コマあたりのミリ秒
t = np.linspace(start_time, end_time, steps)
t_pre = t[t < event_time]
t_post = t[t >= event_time]

strain = np.where(t - event_time >= 0, strain_i, 0)

# ODEの解析
e0 = strain_i       # ODEの引数として入れるためにこの形で定義
s0 = E * strain_i   # バネ要素の応力はステップ歪みに対して即座に応答するため、初期条件として定義
sol = odeint(Maxwell_stepStrain, s0, t_post, args=(tau,)) # ODEの解
stress_pre = np.zeros_like(t_pre)   # ステップ前の応力はゼロ
stress_post = sol[:, 0]             # 応力履歴
stress = np.concatenate([stress_pre, stress_post])

#dt = t[1] - t[0]   # 時間刻み
#integral_stress = np.array([stress[:k+1].sum()*dt for k in range(len(stress))])     # 簡易的な応力の積分
integral_stress = integrate.cumulative_trapezoid(stress, t, initial=0)               # scipyを使った応力の積分

# 描画のためのスケーリング
e = strain/1.0     # 描画のためのスケーリング
s = stress/10**6   # 描画のためのスケーリング ([MPa]単位に変換)
e_s = s/(E/10**6)  # バネの歪み
e_d = integral_stress/eta   # ダッシュポットの歪み
l = 0.1            # [m] 自然長（の半分）
w = 0.5            # 要素の長さと幅の比率
# 直列なので自然長l0=2*lとなっている
el = e*2*l         # [m] 全体の伸び
el_s = e_s*2*l     # [m] バネの伸び
el_d = e_d*2*l     # [m] ダッシュポットの伸び

# グラフの初期設定
fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111, xlabel='$t$ /s')
ax.grid()
title_text = "Maxwell model: step strain"
ax.set_title(title_text)
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
ax.set_xlim(-0.05,0.4)
ax.set_ylim(-5,5)

# 枠組み描画の準備
y_0 = [0, 0]
ax.plot([0, 0.08*l],y_0, c='b')
ax.plot(0,0,'ro', markersize='10')
ax.plot([0,2*l],[-2,-2], c='g')
ax.plot([0,0],[-1.8,-2.2], c='g')
ax.plot([2*l,2*l],[-1.8,-2.2], c='g')

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
var_text = r'$\epsilon_0$ = {0:.1f}, $E$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s'.format(strain_i,E/10**6,eta/10**3)
ax.text(0.5, 0.9, var_text, transform=ax.transAxes)
eq_text = r'd$\sigma$/d$t$ = -$\sigma$/$\tau$'
ax.text(0.5, 0.8, eq_text, transform=ax.transAxes)
res_text = r'$\tau$ = {0:.1f} s'.format(tau)
ax.text(0.5, 0.7, res_text, transform=ax.transAxes)
ax.text(0.35, 0.25, '$l_0$', transform=ax.transAxes)

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

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

# アニメーション更新関数
def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return rod, point, triangle, rod_da, damper, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    # 枠組み
    x_o = 7*l/4 + el[i]
    x_rod = [x_o, l/4+x_o]
    rod.set_data(x_rod,y_0)
    point.set_data([l/4+x_o],[0])
    # ダッシュポット
    x_da = l/2 + el_d[i]
    x_rod_da = [x_da, x_da+3*l/4]
    x_damp = x_da
    y_damp = 0.7*w
    x_damper = [x_damp, x_damp]
    y_damper = [y_damp, -y_damp]
    rod_da.set_data(x_rod_da,y_0)
    damper.set_data(x_damper,y_damper)
    time_text.set_text(time_template % (t[i]))
    # バネ
    x_sp = x_o-el_s[i]-l/2
    x_tri = np.linspace(x_sp, x_o,100)
    y_tri = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri - x_sp)/(x_o-x_sp)-np.pi/2+0.1))-1)
    triangle.set_data(x_tri,y_tri)
    return rod, point, triangle, rod_da, damper, time_text

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

savefile = "./mp4/Maxwell_stepStrain_ani.mp4"
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '30'])

plt.show()