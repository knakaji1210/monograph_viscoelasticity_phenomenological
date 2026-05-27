# 拡張一般化マクスウェルモデルの常微分方程式（ステップ歪み）（アニメーション付き）

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib.animation as animation

'''
テキストの式(6.17)をベースに組み立てる
'''

# 変数の設定
def reqParams():
    E_list = []
    eta_list = []
    tau_list = []
    try:
        infMod = float(input('Enter equilibrium modulus value (MPa) (default = 0.5 MPa): '))*10**6
    except ValueError:
        infMod = 5*10**5
    numComp = 3
    for i in range(numComp):
        try:
            Ej = float(input('Enter modulus value of Maxwell component (MPa) (default = 1 MPa): '))*10**6
        except ValueError:
            Ej = 10**6
        E_list.append(Ej)
        try:
            etaj = float(input('Enter viscosity value of Maxwell component (kPa s) (default = 1000 kPa s): '))*10**3
        except ValueError:
            etaj = 10**6
        eta_list.append(etaj)
        tauj = etaj/Ej
        tau_list.append(tauj)
    return numComp, infMod, E_list, eta_list, tau_list

numComp, infMod, E_list, eta_list, tau_list = reqParams()
insMod = sum(E_list) + infMod  # [Pa] 瞬間弾性率

# 初期条件の設定
try:
    strain_i = float(input('step strain (default = 0.2): '))
except ValueError:
    strain_i = 0.2         # [Pa] ステップ歪み

# ODE解析で用いる関数の定義
def genMaxwell_stepStrain(s, t, e, tau):
# e: 歪み, s: 応力, tau: 緩和時間
# ここでは下でargsとしてe=strain_iを入れてステップ歪みを実現
    dsdt =  - s/tau       # (6.17)
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
stress_array = np.zeros((numComp,len(t_pre)+len(t_post)))  # 各マクスウェル要素の応力履歴を格納する配列
e0 = strain_i
s0 = insMod*e0      # ステップ歪みを加えた直後に応力は全てのバネ要素かかるため、初期条件として定義
for j in range(numComp):    # 各マクスウェル要素の応力を計算
    sj = E_list[j]*e0       # ODEの引数として入れるためにこの形で定義
    sol = odeint(genMaxwell_stepStrain, sj, t_post, args=(e0, tau_list[j]))
    stress_pre = np.zeros_like(t_pre)   # ステップ前の応力はゼロ
    stress_post = sol[:, 0]             # 応力履歴
    stress = np.concatenate([stress_pre, stress_post])
    stress_array[j] = stress
stress_total = np.sum(stress_array, axis=0) + infMod*strain  # 各マクスウェル要素の応力と単独バネの応力を合計

# 描画のためのスケーリング
e = strain/1.0                  # 描画のためのスケーリング
s = stress_total/10**6          # 描画のためのスケーリング ([MPa]単位に変換)
s1 = stress_array[0]/10**6      # 各マクスウェル要素の応力履歴（描画のためのスケーリング）
s2 = stress_array[1]/10**6
s3 = stress_array[2]/10**6
ss = infMod*strain/10**6        # 単独バネの応力履歴（描画のためのスケーリング）
e1 = stress_array[0]/E_list[0]  # 各マクスウェル要素のバネの歪み履歴（描画のためのスケーリング）
e2 = stress_array[1]/E_list[1]
e3 = stress_array[2]/E_list[2]
l = 0.1            # [m] 自然長
w = 0.5            # 要素の長さと幅の比率
# 直列なので自然長l0=2*lとなっている
el = e*2*l         # [m] 全体の伸び
el_s1 = e1*2*l     # [m] Maxwell要素#1のバネの伸び
el_d1 = (e-e1)*2*l # [m] Maxwell要素#1のダッシュポットの伸び
el_s2 = e2*2*l     # [m] Maxwell要素#2のバネの伸び
el_d2 = (e-e2)*2*l # [m] Maxwell要素#2のダッシュポットの伸び
el_s3 = e3*2*l     # [m] Maxwell要素#3のバネの伸び
el_d3 = (e-e3)*2*l # [m] Maxwell要素#3のダッシュポットの伸び

# グラフの初期設定
fig = plt.figure(figsize=(8,7), tight_layout=True)
ax = fig.add_subplot(111)
ax.grid()
title_text = "generalized Maxwell model: step strain"
ax.set_title(title_text)
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
ax.set_xlim(-0.05,0.3)
ax.set_ylim(-5,5)

# 枠組みの描画
y_0 = [0, 0]
y_1 = [3, 3]        # マクスウェル要素#1の位置
y_2 = [1, 1]        # マクスウェル要素#2の位置
y_3 = [-1, -1]      # マクスウェル要素#3の位置
y_4 = [-3, -3]      # 単独バネの位置
ax.plot([-l/4,0],y_0, c='b')
ax.plot([0,0],[3,-3], c='b')
ax.plot([0, 0.08*l],y_1, c='b')
ax.plot([0, 0.08*l],y_2, c='b')
ax.plot([0, 0.08*l],y_3, c='b')
ax.plot([0, 3*l/4],y_4, c='b')            # 7*l/8ではないかも
ax.plot([0,2*l],[-4,-4], c='g')
ax.plot([0,0],[-3.8,-4.2], c='g')
ax.plot([2*l,2*l],[-3.8,-4.2], c='g')
ax.plot(-l/4,0,'ro', markersize='10')

# ダッシュポット描画の準備
for j in range(numComp):
    x_d1 = [0.08*l, 0.92*l]
    y_d1 = [w+3-2*j,w+3-2*j]
    y_d2 = [-w+3-2*j, -w+3-2*j]
    ax.plot(x_d1,y_d1, c='b')
    ax.plot(x_d1,y_d2, c='b')
    ax.plot([0.08*l,0.08*l],[w+3-2*j,-w+3-2*j], c='b')
    rect = patches.Rectangle(xy=(0.08*l, -w+3-2*j), width=0.83*l, height=2*w, facecolor='y')
    ax.add_patch(rect)

# テキスト描画
eq_text = r'd$\sigma_j$/d$t$ = - $\sigma_j$/$\tau_j$'
ax.text(0.25, 0.9, eq_text, transform=ax.transAxes)
res_text = r'$\epsilon_0$ = {0:.2f}, $E_i$ = {1:.2f} MPa, $E_\infty$ = {2:.2f} MPa'.format(e0,insMod/10**6, infMod/10**6)
ax.text(0.5, 0.9, res_text, transform=ax.transAxes)
ax.text(0.4, 0.05, '$l_0$', transform=ax.transAxes)

# 枠組みの描画2
bar, = ax.plot([],[], 'b', animated=True)
rod0, = ax.plot([],[], 'b', animated=True)
rod1, = ax.plot([],[], 'b', animated=True)
rod2, = ax.plot([],[], 'b', animated=True)
rod3, = ax.plot([],[], 'b', animated=True)
rod4, = ax.plot([],[], 'b', animated=True)
point, = ax.plot([],[], 'ro', markersize='10', animated=True)

# マクスウェル要素のバネ
triangle1, = ax.plot([],[], 'b', animated=True)
triangle2, = ax.plot([],[], 'b', animated=True)
triangle3, = ax.plot([],[], 'b', animated=True)
# マクスウェル要素のダッシュポット
rod_da1, = ax.plot([],[], 'b', animated=True)
rod_da2, = ax.plot([],[], 'b', animated=True)
rod_da3, = ax.plot([],[], 'b', animated=True)
damper1, = ax.plot([],[], 'b', lw=4, animated=True)
damper2, = ax.plot([],[], 'b', lw=4, animated=True)
damper3, = ax.plot([],[], 'b', lw=4, animated=True)
# 単独のバネ
triangle4, = ax.plot([],[], 'b', animated=True)
# ここでは[],[]としているが、下で***.set_dataで実際の値を入れている

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

# アニメーション更新関数
def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return bar, rod0, rod1, rod2, rod3, rod4, point, triangle1, triangle2, triangle3, rod_da1, rod_da2, rod_da3, damper1, damper2, damper3, triangle4, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    # 枠組み
    x_o = 2*l + el[i]
    bar.set_data([x_o,x_o],[3,-3])
    x_rod0 = [x_o, l/4+x_o]
    x_rod1 = [-l/4+x_o, x_o]
    x_rod2 = [-l/4+x_o, x_o]
    x_rod3 = [-l/4+x_o, x_o]
    x_rod4 = [-3*l/4+x_o, x_o]
    rod0.set_data(x_rod0,y_0)
    rod1.set_data(x_rod1,y_1)
    rod2.set_data(x_rod2,y_2)
    rod3.set_data(x_rod3,y_3)
    rod4.set_data(x_rod4,y_4)
    point.set_data([l/4+x_o],[0])

    # マクスウェル要素#1のダッシュポット
    x_da1 = l/2 + el_d1[i]
    x_rod_da1 = [x_da1, x_da1+3*l/4]
    x_damp1 = x_da1
    y_damp = 0.7*w
    x_damper1 = [x_damp1, x_damp1]
    y_damper1 = [y_damp+3, -y_damp+3]
    rod_da1.set_data(x_rod_da1,y_1)
    damper1.set_data(x_damper1,y_damper1)

    # マクスウェル要素要素#1のバネ
    x_sp1 = x_da1+3*l/4
    x_tri1 = np.linspace(x_sp1, x_o-l/4,100)
    y_tri1 = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri1 - x_sp1)/(x_o - l/4 -x_sp1)-np.pi/2+0.1))+5)
    triangle1.set_data(x_tri1,y_tri1)

    # マクスウェル要素#2のダッシュポット
    x_da2 = l/2 + el_d2[i]
    x_rod_da2 = [x_da2, x_da2+3*l/4]
    x_damp2 = x_da2
    y_damp = 0.7*w
    x_damper2 = [x_damp2, x_damp2]
    y_damper2 = [y_damp+1, -y_damp+1]
    rod_da2.set_data(x_rod_da2,y_2)
    damper2.set_data(x_damper2,y_damper2)

    # マクスウェル要素要素#2のバネ
    x_sp2 = x_da2+3*l/4
    x_tri2 = np.linspace(x_sp2, x_o-l/4,100)
    y_tri2 = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri2 - x_sp2)/(x_o - l/4 -x_sp2)-np.pi/2+0.1))+1)
    triangle2.set_data(x_tri2,y_tri2)

    # マクスウェル要素#3のダッシュポット
    x_da3 = l/2 + el_d3[i]
    x_rod_da3 = [x_da3, x_da3+3*l/4]
    x_damp3 = x_da3
    y_damp = 0.7*w
    x_damper3 = [x_damp3, x_damp3]
    y_damper3 = [y_damp-1, -y_damp-1]
    rod_da3.set_data(x_rod_da3,y_3)
    damper3.set_data(x_damper3,y_damper3)

    # マクスウェル要素要素#3のバネ
    x_sp3 = x_da3+3*l/4
    x_tri3 = np.linspace(x_sp3, x_o-l/4,100)
    y_tri3 = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri3 - x_sp3)/(x_o - l/4 -x_sp3)-np.pi/2+0.1))-3)
    triangle3.set_data(x_tri3,y_tri3)

    # 単独のバネ
    x_sp2 = -3*l/4+x_o
    x_tri2 = np.linspace(3*l/4, x_sp2,100)
    y_tri2 = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri2 - 3*l/4)/(x_sp2-3*l/4)-np.pi/2+0.1))-7)
    triangle4.set_data(x_tri2,y_tri2)
    time_text.set_text(time_template % (t[i]))
    return bar, rod0, rod1, rod2, rod3, rod4, point, triangle1, triangle2, triangle3, rod_da1, rod_da2, rod_da3, damper1, damper2, damper3, triangle4, time_text

'''
y_triの中の重要部分は
x_tri1 = np.linspace(a, b,100)
のとき
(xtri - a)/(b - a)
になる 
'''

# テーブルの追加（グラフの下に配置）
param_table = ax.table(
    cellText=[['{0:.1f}'.format(np.log10(E_list[j])) for j in range(numComp)],
              ['{0:.1f}'.format(eta_list[j]/10**3) for j in range(numComp)],
              ['{0:.2f}'.format(np.log10(tau_list[j])) for j in range(numComp)]],
    rowLabels=[r"log($E_j$ /Pa)", r"$\eta_j$ /kPa s", r"log($\tau_j$ /s)"],
    colLabels=['#{}'.format(j+1) for j in range(numComp)],
    loc='bottom',
    bbox=[0.0, -0.45, 1.0, 0.3]  # [x, y, 幅, 高さ] で微調整
)
param_table.set_fontsize(11)
param_table.scale(1, 1.5) # セルの大きさを調整
plt.subplots_adjust(bottom=0.2)

# アニメーション実行   
ani = animation.FuncAnimation(fig, update, frames=steps, 
                    init_func=init, blit=True, interval=interval_ms, repeat=False)

savefile = './mp4/genMaxwell_stepStrain_ani.mp4'
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '30'])

plt.tight_layout()
plt.show()