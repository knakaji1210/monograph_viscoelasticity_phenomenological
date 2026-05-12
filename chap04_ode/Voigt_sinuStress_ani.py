# フォークトモデルの常微分方程式（振動応力）（アニメーション付き）

import numpy as np
from scipy.integrate import odeint
from scipy import integrate
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib.animation as animation

'''
テキストの式(2.23)をベースに組み立てる
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

tau = eta/E                     # [s] 遅延時間

# 振動応力の設定
try:
    samp = float(input('amplitude for sinusoidal stress [MPa] (default=0.04): '))*10**6
except ValueError:
    samp = 0.04*10**6
try:
    freq = float(input('frequency for sinusoidal stress [Hz] (default=0.3813 Hz): '))
except ValueError:
    freq = 1/np.pi

af = 2*np.pi*freq

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
start_time = np.where(1.0/freq > 2.0, -2.0, -1.0/freq)      # 開始時間
#print('start_time = {0:.2f} s'.format(start_time))  
end_time = np.where(4.0/freq > 20.0, 2.0/freq, 4.0/freq)    # 終了時間
#print('end_time = {0:.2f} s'.format(end_time))
event_time = 0.0    # ステップ応力を加える時刻
time_duration = end_time - start_time       # [s] 継続時間
time_duration_pre = event_time - start_time # [s] ステップ前の継続時間
time_duration_post = end_time - event_time  # [s] ステップ後の継続時間
fps = 60            # 1秒あたりのフレーム数、30だと足りないので60に変更
steps = int(time_duration * fps) + 1        # 総フレーム数
interval_ms = 1000 / fps                    # 1コマあたりのミリ秒
t = np.linspace(start_time, end_time, steps)
t_pre = t[t < event_time]
t_post = t[t >= event_time]

stress_post = samp*np.sin(af*t)  # 振動応力の関数
stress = np.where(t - event_time >= 0, stress_post, 0)

# ODEの解析
e0 = 0                              # 初期条件として定義
sol = odeint(Voigt_sinuStress, e0, t_post, args=(samp,af,E,tau)) # ODEの解
strain_pre = np.zeros_like(t_pre)   # ステップ前の歪みはゼロ
strain_post = sol[:, 0]             # 歪み履歴
strain = np.concatenate([strain_pre, strain_post])

dedt = np.gradient(strain, t) 

# 位相差の計算
stress_latter = stress[int(0.4*len(stress)):]     # 後半部分を抽出（前半は過渡応答を含むから）
strain_latter = strain[int(0.4*len(strain)):]     # 後半部分を抽出（前半は過渡応答を含むから）
strain_max = np.max(strain_latter)
ind = getNearestIndex2value(strain_latter,0)      # 出力信号が0になるindexを抽出           
phase_diff = (180/np.pi)*np.arcsin(np.abs(stress_latter[ind])/samp)

# 描画のためのスケーリング
s = stress/10**6        # 描画のためのスケーリング ([MPa]単位に変換)
e = strain/1.0          # 描画のためのスケーリング
e_max = strain_max/1.0
s_s = E*e/10**6         # バネの応力 ([MPa]単位に変換)
s_d = eta*dedt/10**6    # ダッシュポットの応力 ([MPa]単位に変換)
l = 0.1                 # [m] 自然長
w = 0.5                 # 要素の長さと幅の比率
# 並列なので自然長l0=lとなっている
el = e*l                # [m] 全体の伸び

# グラフの初期設定
fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111, xlabel='$t$ /s')
ax.grid()
title_text = "Voigt model: sinusoidal stress"
ax.set_title(title_text)
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
ax.set_xlim(-0.05,0.3)
ax.set_ylim(-5,5)

# 枠組み描画の準備
y_0 = [0, 0]
y_1 = [1, 1]
y_2 = [-1, -1]
ax.plot([-l/4, 0],y_0, c='b')
ax.plot([0, l],[-4,-4], c='g')
ax.plot([0, 0],[-3.8,-4.2], c='g')
ax.plot([l, l],[-3.8,-4.2], c='g')
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
var_text = r'$\sigma_{{amp}}$ = {0:.2f} MPa, $f$ = {1:.3f} Hz, $E$ = {2:.1f} MPa, $\eta$ = {3:.1f} kPa s'.format(samp/10**6,freq,E/10**6,eta/10**3)
ax.text(0.4, 0.9, var_text, transform=ax.transAxes)
eq_text = r'd$\epsilon$/d$t$ = ($\sigma$/$E$ - $\epsilon$)/$\tau$'
ax.text(0.4, 0.8, eq_text, transform=ax.transAxes)
res_text = r'$\tau$ = {0:.1f} s'.format(tau)
ax.text(0.4, 0.7, res_text, transform=ax.transAxes)
ax.text(0.30, 0.05, '$l_0$', transform=ax.transAxes)
ax.text(0.75, 0.28, '$\sigma$ (input)', transform=ax.transAxes)
ax.text(0.75, 0.18, '$\epsilon$ (output)', transform=ax.transAxes)
af_text = r'$\omega$ = {0:.2f} s$^{{-1}}$'.format(af)
ax.text(0.1, 0.8, af_text, transform=ax.transAxes)
aft_text = r'$\omega\tau$ = {0:.2f}'.format(af*tau)
ax.text(0.1, 0.7, aft_text, transform=ax.transAxes)
eamp_template = r'$\epsilon_{{amp}}$ = {0:.2f}'.format(e_max)
ax.text(0.75, 0.42, eamp_template, transform=ax.transAxes)
phase_diff_template = r'$\theta$ = {0:.1f} $\degree$'.format(phase_diff)
ax.text(0.75, 0.35, phase_diff_template, transform=ax.transAxes)

# 枠組みの描画
bar, = ax.plot([],[], 'b', animated=True)
rod, = ax.plot([],[], 'b', animated=True)
point, = ax.plot([],[], 'ro', markersize='10', animated=True)

# バネの描画
rod_sp, = ax.plot([],[], 'b', animated=True)
triangle, = ax.plot([],[], 'b', animated=True)

# ダッシュポットの描画
rod_da, = ax.plot([],[], 'b', animated=True)
damper, = ax.plot([],[], 'b', lw=4, animated=True)
# ここでは[],[]としているが、下で***.set_dataで実際の値を入れている

# 応力の描画
stress, = ax.plot([],[], 'r', lw=2, animated=True)
arrow_s_p, = ax.plot([],[], 'r', marker=9, markersize='10', animated=True)
arrow_s_n, = ax.plot([],[], 'r', marker=8, markersize='10', animated=True)

# 歪みの描画
strain, = ax.plot([],[], 'b', lw=2, animated=True)
arrow_e_p, = ax.plot([],[], 'b', marker=9, markersize='10', animated=True)
arrow_e_n, = ax.plot([],[], 'b', marker=8, markersize='10', animated=True)

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

# アニメーション更新関数
def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return bar, rod, point, rod_sp, triangle, rod_da, damper, stress, arrow_s_p, arrow_s_n, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    # 枠組み
    x_o = 3*l/4 + el[i]
    x_c = l/4 + x_o
    bar.set_data([x_c,x_c],[1,-1])
    x_rod = [x_c, l/4+x_c]
    rod.set_data(x_rod,y_0)
    point.set_data([l/4+x_c],[0])
    # バネ
    x_rod_sp = [x_o, x_c]
    rod_sp.set_data(x_rod_sp,y_2)
    x_tri = np.linspace(l/4, x_o,100)
    y_tri = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri - l/4)/(el[i]+l/2)-np.pi/2+0.1))-3)
    triangle.set_data(x_tri,y_tri)
    # ダッシュポット
    x_rod_da = [x_o-l/4, x_c]
    x_damp = x_o-l/4
    y_damp = 0.7*w
    x_damper = [x_damp, x_damp]
    y_damper = [y_damp+1, -y_damp+1]
    rod_da.set_data(x_rod_da,y_1)
    damper.set_data(x_damper,y_damper)
    # 応力
    a = 0.2    # 見かけ上の振幅
    x_stress = [5*l/4, 5*l/4 + a*s[i]]
    stress.set_data(x_stress,[-2,-2])
    if s[i] > 0:
        arrow_s_p.set_data([5*l/4 + a*s[i]],[-2])
    else:
        arrow_s_n.set_data([5*l/4 + a*s[i]],[-2])
    # 歪み
    x_strain = [5*l/4, 5*l/4 + el[i]]
    strain.set_data(x_strain,[-3,-3])
    if el[i] > 0:
        arrow_e_p.set_data([5*l/4 + el[i]],[-3])
    else:
        arrow_e_n.set_data([5*l/4 + el[i]],[-3])
    time_text.set_text(time_template % (t[i]))
    return bar, rod, point, rod_sp, triangle, rod_da, damper, stress, arrow_s_p, arrow_s_n, time_text

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

savefile = './mp4/Voigt_sinuStress_ani_(f={0:.2f}Hz).mp4'.format(freq)
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '60'])

plt.tight_layout()
plt.show()