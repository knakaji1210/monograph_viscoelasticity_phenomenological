# SLS IIモデルの常微分方程式（ステップ歪み）（アニメーション付き）

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib.animation as animation

'''
テキストの式(5.40)をベースに組み立てる
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
    eta = float(input('viscosity [kPa s] (default = 500.0 kPa s): '))*10**3
except ValueError:
    eta = 5*10**5               # [Pa s] 粘度

insMod = E1+E2              # [Pa] 瞬間弾性率
infMod = E2                 # [Pa] 緩和弾性率
k = insMod/infMod
tau = eta/E1                # [s] 緩和時間

# 初期条件の設定
try:
    strain_i = float(input('step strain (default = 0.2): '))
except ValueError:
    strain_i = 0.2         # [Pa] ステップ歪み

# ODE解析で用いる関数の定義
def SLS2_stepStrain(s, t, e, infMod, tau):
# e: 歪み, s: 応力, infMod: 緩和弾性率, tau: 緩和時間
# ここでは下でargsとしてe=strain_iを入れてステップ歪みを実現
    dsdt = (infMod*e - s)/tau   # (5.40)
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
e0 = strain_i   # ODEの引数として入れるためにこの形で定義
s0 = insMod*e0      # # ステップ歪みを加えた直後に応力はspring1とspring2の両方にかかるため、初期条件として定義
sol = odeint(SLS2_stepStrain, s0, t_post, args=(e0,infMod,tau))
stress_pre = np.zeros_like(t_pre)   # ステップ前の応力はゼロ
stress_post = sol[:, 0]             # 応力履歴
stress = np.concatenate([stress_pre, stress_post])

# 描画のためのスケーリング
e = strain/1.0          # 描画のためのスケーリング
s = stress/10**6        # 描画のためのスケーリング ([MPa]単位に変換)
s2 = E2*e/10**6         # [MPa] バネ2の応力
s1 = s - s2             # [MPa] マクスウェル要素の応力
e1 = s1/E1*10**6        # [] バネ1の歪み
e2 = e - e1             # [] ダッシュポットの歪み
de2dt = np.gradient(e2, t)  # numpyを使ったダッシュポットの歪みの微分
s1_d = eta*de2dt/10**6  # [MPa] ダッシュポットの応力（確認用:s1と一致するはず）
l = 0.1            # [m] 自然長
w = 0.5            # 要素の長さと幅の比率
# 直列なので自然長l0=2*lとなっている
el = e*2*l         # [m] 全体の伸び
el_s = e1*2*l      # [m] Maxwell要素のバネ1の伸び
el_d = e2*2*l      # [m] Maxwell要素のダッシュポットの伸び

# グラフの初期設定
fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111)
ax.grid()
title_text = "SLS II model: step strain"
ax.set_title(title_text)
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
ax.set_xlim(-0.05,0.3)
ax.set_ylim(-5,5)

# 枠組みの描画
y_0 = [0, 0]
y_1 = [1, 1]
y_2 = [-1, -1]
ax.plot([-l/4,0],y_0, c='b')
ax.plot([0,0],[1,-1], c='b')
ax.plot([0, 0.08*l],y_1, c='b')
ax.plot([0, 3*l/4],y_2, c='b')            # 7*l/8ではないかも
ax.plot([0,2*l],[-3,-3], c='g')
ax.plot([0,0],[-2.8,-3.2], c='g')
ax.plot([2*l,2*l],[-2.8,-3.2], c='g')
ax.plot(-l/4,0,'ro', markersize='10')

# ダッシュポット描画の準備
x_d1 = [0.08*l, 0.92*l]
y_d1 = [w+1, w+1]
y_d2 = [-w+1, -w+1]
ax.plot(x_d1,y_d1, c='b')
ax.plot(x_d1,y_d2, c='b')
ax.plot([0.08*l,0.08*l],[w+1,-w+1], c='b')
rect = patches.Rectangle(xy=(0.08*l, -w+1), width=0.83*l, height=2*w, facecolor='y')
ax.add_patch(rect)

# テキスト描画
var_text = r'$\epsilon_0$ = {0:.2f}, $E_1$ = {1:.1f} MPa, $E_2$ = {2:.1f} MPa, $\eta$ = {3:.1f} kPa s'.format(e0,E1/10**6,E2/10**6,eta/10**3)
ax.text(0.4, 0.9, var_text, transform=ax.transAxes)
eq_text = r'd$\sigma$/d$t$ = ($E_{{{\infty}}}$$\epsilon_0$ - $\sigma$)/$\tau$'
ax.text(0.4, 0.8, eq_text, transform=ax.transAxes)
res_text = r'$E_i$ = {0:.2f} MPa, $E_\infty$ = {1:.2f} MPa, $\tau$ = {2:.2f} s'.format(insMod/10**6, infMod/10**6, tau)
ax.text(0.4, 0.7, res_text, transform=ax.transAxes)
ax.text(0.4, 0.15, '$l_0$', transform=ax.transAxes)

# 枠組みの描画2
bar, = ax.plot([],[], 'b', animated=True)
rod1, = ax.plot([],[], 'b', animated=True)
rod2, = ax.plot([],[], 'b', animated=True)
rod3, = ax.plot([],[], 'b', animated=True)
point, = ax.plot([],[], 'ro', markersize='10', animated=True)

# マクスウェル要素のバネ
triangle1, = ax.plot([],[], 'b', animated=True)
# マクスウェル要素のダッシュポット
rod_da, = ax.plot([],[], 'b', animated=True)
damper, = ax.plot([],[], 'b', lw=4, animated=True)
# 単独のバネ
triangle2, = ax.plot([],[], 'b', animated=True)
# ここでは[],[]としているが、下で***.set_dataで実際の値を入れている

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

# アニメーション更新関数
def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return bar, rod1, rod2, rod3, point, triangle1, rod_da, damper, triangle2, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    # 枠組み
    x_o = 2*l + el[i]
    bar.set_data([x_o,x_o],[1,-1])
    x_rod1 = [x_o, l/4+x_o]
    x_rod2 = [-l/4+x_o, x_o]
    x_rod3 = [-3*l/4+x_o, x_o]
    rod1.set_data(x_rod1,y_0)
    rod2.set_data(x_rod2,y_1)
    rod3.set_data(x_rod3,y_2)
    point.set_data([l/4+x_o],[0])

    # マクスウェル要素要素のダッシュポット
    x_da = l/2 + el_d[i]
    x_rod_da = [x_da, x_da+3*l/4]
    x_damp = x_da
    y_damp = 0.7*w
    x_damper = [x_damp, x_damp]
    y_damper = [y_damp+1, -y_damp+1]
    rod_da.set_data(x_rod_da,y_1)
    damper.set_data(x_damper,y_damper)

    # マクスウェル要素要素のバネ
    x_sp1 = x_da+3*l/4
    x_tri1 = np.linspace(x_sp1, x_o-l/4,100)
    y_tri1 = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri1 - x_sp1)/(x_o - l/4 -x_sp1)-np.pi/2+0.1))+1)
    triangle1.set_data(x_tri1,y_tri1)

    # 単独のバネ
    x_sp2 = -3*l/4+x_o
    x_tri2 = np.linspace(3*l/4, x_sp2,100)
    y_tri2 = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri2 - 3*l/4)/(x_sp2-3*l/4)-np.pi/2+0.1))-3)
    triangle2.set_data(x_tri2,y_tri2)
    time_text.set_text(time_template % (t[i]))
    return bar, rod1, rod2, rod3, point, triangle1, rod_da, damper, triangle2, time_text

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

savefile = './mp4/SLS2_stepStrain_ani.mp4'
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '30'])

plt.tight_layout()
plt.show()