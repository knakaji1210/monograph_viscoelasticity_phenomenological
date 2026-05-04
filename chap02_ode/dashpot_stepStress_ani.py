# ダッシュポット要素の常微分方程式（ステップ応力）（アニメーション付き）

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib.animation as animation

'''
テキストの式(2.6)をベースに組み立てる
'''

# 変数の設定
try:
    eta = float(input('viscosity [kPa s] (default = 500.0 kPa s): '))*10**3
except ValueError:
    eta = 5*10**5           # [Pa s] 粘度

# 初期条件の設定
try:
    stress_i = float(input('step stress [MPa] (default = 0.02 MPa): '))*10**6
except ValueError:
    stress_i = 0.02*10**6   # [Pa] ステップ応力

# ODE解析で用いる関数の定義
def dashpot_stepStress(e, t, s, eta):
# e: 歪み, s: 応力, eta: 粘度
# ここでは下でargsとしてs0=stress_iを入れてステップ応力を実現
    dedt = s/eta    # (2.6)
    return dedt

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

stress = np.where(t - event_time >= 0, stress_i, 0)

# ODEの解析
s0 = stress_i   # ODEの引数として入れるためにこの形で定義
e0 = 0.0        # ステップ応力を加える前の歪みはゼロとするため、初期条件として定義
sol = odeint(dashpot_stepStress, e0, t_post, args=(s0,eta)) # ODEの解
strain_pre = np.zeros_like(t_pre)   # ステップ前の歪みはゼロ
strain_post = sol[:, 0]             # 歪み履歴
strain = np.concatenate([strain_pre, strain_post])

# 描画のためのスケーリング
e = strain/1.0     # 描画のためのスケーリング
s = stress/10**6   # 描画のためのスケーリング ([MPa]単位に変換)
l = 0.1            # [m] 自然長
w = 0.5            # 要素の長さと幅の比率
el = e*l           # [m] 全体の伸び

# グラフの初期設定
fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111)
ax.grid()
title_text = "dashpot (Newton's viscosity): step stress"
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
ax.plot([0,l],[-2,-2], c='g')
ax.plot([0,0],[-1.8,-2.2], c='g')
ax.plot([l,l],[-1.8,-2.2], c='g')

# テキスト描画
var_text = r'$\sigma_0$ = {0:.2f} MPa, $\eta$ = {1:.1f} kPa s'.format(s0/10**6,eta/10**3)
ax.text(0.5, 0.9, var_text, transform=ax.transAxes)
eq_text = r'd$\epsilon$/d$t$ = $\sigma_0$/$\eta$'
ax.text(0.5, 0.8, eq_text, transform=ax.transAxes)
ax.text(0.3, 0.25, '$l_0$', transform=ax.transAxes)

# ダッシュポットの描画
rod, = ax.plot([],[], 'b', animated=True)
damper, = ax.plot([],[], 'b', lw=4, animated=True)
point, = ax.plot([], [], 'ro', markersize='10', animated=True)
# ここでは[],[]としているが、下で***.set_dataで実際の値を入れている

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

# アニメーション更新関数
def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return rod, damper, point, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    x_rod = [l/2 + el[i], l + el[i]]
    x_damp = l/2+el[i]
    y_damp = 0.7*w
    x_damper = [x_damp, x_damp]
    y_damper = [y_damp, -y_damp]
    rod.set_data(x_rod,y_0)
    damper.set_data(x_damper,y_damper)
    point.set_data([l + el[i]],[0])
    time_text.set_text(time_template % (t[i]))
    return rod, damper, point, time_text

# アニメーション実行
ani = animation.FuncAnimation(fig, update, frames=steps, 
                    init_func=init, blit=True, interval=interval_ms, repeat=False)

savefile = "./mp4/dashpot_stepStress_ani.mp4"
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '30'])

plt.show()