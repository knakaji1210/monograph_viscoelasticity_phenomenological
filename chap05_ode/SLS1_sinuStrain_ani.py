# SLS Iモデルの常微分方程式（振動歪み）（アニメーション付き）

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib.animation as animation

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
    eamp = float(input('amplitude for sinusoidal strain [] (default=0.1): '))
except ValueError:
    eamp = 0.1
try:
    freq = float(input('frequency for sinusoidal strain [Hz] (default=1.91 Hz): '))
except ValueError:
    freq = 6/np.pi

af = 2*np.pi*freq

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
start_time = np.where(1.0/freq > 2.0, -2.0, -1.0/freq)      # 開始時間
#print('start_time = {0:.2f} s'.format(start_time))  
end_time = np.where(4.0/freq > 20.0, 2.0/freq, 4.0/freq)    # 終了時間
#print('end_time = {0:.2f} s'.format(end_time))
event_time = 0.0    # ステップ歪みを加える時刻
time_duration = end_time - start_time       # [s] 継続時間
time_duration_pre = event_time - start_time # [s] ステップ前の継続時間
time_duration_post = end_time - event_time  # [s] ステップ後の継続時間
fps = 60            # 1秒あたりのフレーム数、30だと足りないので60に変更
#fps = 500           # 高周波数の場合
steps = int(time_duration * fps) + 1        # 総フレーム数
interval_ms = 1000 / fps                    # 1コマあたりのミリ秒
t = np.linspace(start_time, end_time, steps)
t_pre = t[t < event_time]
t_post = t[t >= event_time]

strain_post = eamp*np.sin(af*t)  # 振動歪みの関数
strain = np.where(t - event_time >= 0, strain_post, 0)

# ODEの解析
s0 = 0                              # 初期条件として定義
sol = odeint(SLS1_sinuStrain, s0, t_post, args=(eamp,af,insMod,k,tau))
stress_pre = np.zeros_like(t_pre)   # ステップ前の応力はゼロ
stress_post = sol[:, 0]             # 応力履歴
stress = np.concatenate([stress_pre, stress_post])

# 位相差の計算
strain_latter = strain[int(0.4*len(strain)):]     # 後半部分を抽出（前半は過渡応答を含むから）
stress_latter = stress[int(0.4*len(stress)):]     # 後半部分を抽出（前半は過渡応答を含むから）
stress_max = np.max(stress_latter)
ind = getNearestIndex2value(stress_latter,0)      # 出力信号が0になるindexを抽出           
phase_diff = (180/np.pi)*np.arcsin(np.abs(strain_latter[ind])/eamp)

# 描画のためのスケーリング
e = strain/1.0          # 描画のためのスケーリング
s = stress/10**6        # 描画のためのスケーリング ([MPa]単位に変換)
s_max = stress_max/10**6
e1 = s*10**6/E1         # [] バネ1の歪み
e2 = e - e1             # [] フォークト要素の歪み
s1 = E2*e2/10**6        # [MPa] フォークト要素のバネ2の応力
de2dt = np.gradient(e2, t)  # numpyを使ったフォークト要素の歪みの微分
s2 = eta*de2dt/10**6    # フォークト要素のダッシュポットの応力 ([MPa]単位に変換)
l = 0.1            # [m] 自然長
w = 0.5            # 要素の長さと幅の比率
# 直列なので自然長l0=2*lとなっている
el = e*2*l         # [m] 全体の伸び
el_s = e1*2*l        # [m] 単独バネの伸び
el_v = e2*2*l        # [m] Voigt要素部分の伸び

# グラフの初期設定
fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111)
ax.grid()
title_text = "SLS I model: sinusoidal strain"
ax.set_title(title_text)
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
ax.set_xlim(-0.05,0.3)
ax.set_ylim(-5,5)

# 枠組みの描画
y_0 = [0, 0]
y_1 = [1, 1]
y_2 = [-1, -1]
ax.plot([-l/4, 0],y_0, c='b')
ax.plot([0, 2*l],[-4,-4], c='g')
ax.plot([0, 0],[-3.8,-4.2], c='g')
ax.plot([2*l, 2*l],[-3.8,-4.2], c='g')
ax.plot([0, 0],[1,-1], c='b')
ax.plot([0, 0.2*l],y_1, c='b')
ax.plot([0, l/4],y_2, c='b')
ax.plot(-l/4,0,'ro', markersize='10')

# ダッシュポット描画の準備
x_d1 = [0.2*l, 0.8*l]
y_d1 = [w+1, w+1]
y_d2 = [-w+1, -w+1]
ax.plot(x_d1,y_d1, c='b')
ax.plot(x_d1,y_d2, c='b')
ax.plot([0.2*l,0.2*l],[w+1,-w+1], c='b')
rect = patches.Rectangle(xy=(0.2*l, -w+1), width=0.6*l, height=2*w, facecolor='y')
ax.add_patch(rect)

# テキスト描画
var_text = r'$\epsilon_{{amp}}$ = {0:.2f}, $E_1$ = {1:.1f} MPa, $E_2$ = {2:.1f} MPa, $\eta$ = {3:.1f} kPa s'.format(eamp,E1/10**6,E2/10**6,eta/10**3)
ax.text(0.4, 0.9, var_text, transform=ax.transAxes)
eq_text = r'd$\sigma$/d$t$ = ($E_i$$\epsilon$ + $E_i$$\tau$ d$\epsilon$/d$t$ - $k$$\sigma$)/$\tau$'
ax.text(0.4, 0.8, eq_text, transform=ax.transAxes)
res_text = r'$E_i$ = {0:.2f} MPa, $E_\infty$ = {1:.2f} MPa, $\tau$ = {2:.2f} s'.format(insMod/10**6, infMod/10**6, tau)
ax.text(0.4, 0.7, res_text, transform=ax.transAxes)
ax.text(0.4, 0.05, '$l_0$', transform=ax.transAxes)
ax.text(0.85, 0.28, '$\epsilon$ (input)', transform=ax.transAxes)
ax.text(0.85, 0.18, '$\sigma$ (output)', transform=ax.transAxes)
af_text = r'$\omega$ = {0:.2f} s$^{{-1}}$'.format(af)
ax.text(0.1, 0.8, af_text, transform=ax.transAxes)
aft_text = r'$\omega\tau$ = {0:.2f}'.format(af*tau)
ax.text(0.1, 0.7, aft_text, transform=ax.transAxes)
samp_text = r'$\sigma_{{amp}}$ = {0:.3f}'.format(s_max)
ax.text(0.85, 0.42, samp_text, transform=ax.transAxes)
phase_diff_text = r'$\theta$ = {0:.1f} $\degree$'.format(phase_diff)
ax.text(0.85, 0.35, phase_diff_text, transform=ax.transAxes)

# 枠組みの描画2
bar, = ax.plot([],[], 'b', animated=True)
rod, = ax.plot([],[], 'b', animated=True)
point, = ax.plot([],[], 'ro', markersize='10', zorder=10, animated=True)

# フォークト要素のバネ
rod_sp2, = ax.plot([],[], 'b', animated=True)
triangle2, = ax.plot([],[], 'b', animated=True)
# フォークト要素のダッシュポット
rod_da, = ax.plot([],[], 'b', animated=True)
damper, = ax.plot([],[], 'b', lw=4, animated=True)
# 単独のバネ
rod_sp1, = ax.plot([],[], 'b', animated=True)
triangle1, = ax.plot([],[], 'b', animated=True)
# ここでは[],[]としているが、下で***.set_dataで実際の値を入れている

# 歪みの描画
strain_bar, = ax.plot([],[], 'b', lw=2, animated=True)
arrow_e_p, = ax.plot([],[], 'b', marker=9, markersize='10', animated=True)
arrow_e_n, = ax.plot([],[], 'b', marker=8, markersize='10', animated=True)

# 応力の描画
stress_bar, = ax.plot([],[], 'r', lw=2, animated=True)
arrow_s_p, = ax.plot([],[], 'r', marker=9, markersize='10', animated=True)
arrow_s_n, = ax.plot([],[], 'r', marker=8, markersize='10', animated=True)

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

# アニメーション更新関数
def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return bar, rod, point, rod_sp1, rod_sp2, triangle1, triangle2, rod_da, damper, strain_bar, arrow_e_p, arrow_e_n, stress_bar, arrow_s_p, arrow_s_n, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    # 枠組み
    x_o = l + el_v[i]
    x_p = 2*l+el[i]
    bar.set_data([x_o,x_o],[1,-1])
    x_rod = [x_o, l/4 + x_o]
    rod.set_data(x_rod,y_0)
    # フォークト要素のバネ
    x_rod_sp2 = [x_o - l/4, x_o]
    rod_sp2.set_data(x_rod_sp2,y_2)
    x_tri2 = np.linspace(l/4, x_o - l/4,100)
    y_tri2 = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri2 - l/4)/(el_v[i]+l/2)-np.pi/2+0.1))-3)
    triangle2.set_data(x_tri2,y_tri2)
    # フォークト要素のダッシュポット
    x_rod_da = [x_o-l/2, x_o]
    x_damp = x_o - l/2
    y_damp = 0.7*w
    x_damper = [x_damp, x_damp]
    y_damper = [y_damp+1, -y_damp+1]
    rod_da.set_data(x_rod_da,y_1)
    damper.set_data(x_damper,y_damper)
    # 単独のバネ
    x_rod_sp1 = [x_p-l/4 , x_p]
    rod_sp1.set_data(x_rod_sp1,y_0)
    x_tri1 = np.linspace(x_o+l/4, x_p-l/4,100)
    y_tri1 = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri1 - x_o - l/4)/(el_s[i]+l/2)-np.pi/2+0.1))-1)
    triangle1.set_data(x_tri1,y_tri1)
    point.set_data([x_p],[0])    
    # 歪み
    x_strain = [2*l, 2*l + el[i]]
    strain_bar.set_data(x_strain,[-2,-2])
    if el[i] > 0:
        arrow_e_p.set_data([2*l + el[i]],[-2])
    else:
        arrow_e_n.set_data([2*l + el[i]],[-2])
    # 応力
    a = 0.2    # 見かけ上の振幅
    x_stress = [2*l, 2*l + a*s[i]]
    stress_bar.set_data(x_stress,[-3,-3])
    if s[i] > 0:
        arrow_s_p.set_data([2*l + a*s[i]],[-3])
    else:
        arrow_s_n.set_data([2*l + a*s[i]],[-3])
    time_text.set_text(time_template % (t[i]))
    return bar, rod, point, rod_sp1, rod_sp2, triangle1, triangle2, rod_da, damper, strain_bar, arrow_e_p, arrow_e_n, stress_bar, arrow_s_p, arrow_s_n, time_text

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

savefile = './mp4/SLS1_sinuStrain_ani_(f={0:.2f}Hz).mp4'.format(freq)
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '60'])
#ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '500'])  # 高周波数の場合

plt.tight_layout()
plt.show()