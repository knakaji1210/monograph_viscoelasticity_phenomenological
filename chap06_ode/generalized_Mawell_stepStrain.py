# 拡張一般化マクスウェルモデルの常微分方程式（ステップ歪み）

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation

'''
テキストの式(6.20)をベースに組み立てる
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
    dsdt =  - s/tau       # (6.20)
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

# グラフの初期設定
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7.9), sharex=True, tight_layout=True)

# --- 上段：歪み (Input) ---
ax1.set_xlim(start_time, end_time)
ax1.set_ylim(-np.max(e)*0.5, np.max(e)*1.8) # 縦軸を固定
ax1.set_ylabel('Applied strain, $\epsilon$ /')
ax1.set_title("generalized Maxwell model: step strain")
ax1.grid(True, ls='--')

line_strain, = ax1.plot([], [], color='blue', lw=2, label='Step strain (t = {0:.1f} s)'.format(event_time))
line_strain_e1, = ax1.plot([], [], color='cyan', ls="dashed", lw=1, label='$\epsilon$ (spring in Maxwell #1)')
line_strain_e2, = ax1.plot([], [], color='gray', ls="dashed", lw=1, label='$\epsilon$ (spring in Maxwell #2)')
line_strain_e3, = ax1.plot([], [], color='lightgray', ls="dashed", lw=1, label='$\epsilon$ (spring in Maxwell #3)')
ax1.legend(loc='upper right')

# --- 下段：応力 (Response) ---
ax2.set_xlim(start_time, end_time)
ax2.set_ylim(-np.max(s)*0.5, np.max(s)*1.5) # 縦軸を固定
ax2.set_xlabel('$t$ /s')
ax2.set_ylabel('Stress, $\sigma$ /MPa')
ax2.grid(True, ls='--')

line_stress, = ax2.plot([], [], color='red', lw=2, label='Response to step strain')
line_stress_s1, = ax2.plot([], [], color='green', ls="dashed", lw=1, label='$\sigma$ (Maxwell #1)')
line_stress_s2, = ax2.plot([], [], color='orange', ls="dashed", lw=1, label='$\sigma$ (Maxwell #2)')
line_stress_s3, = ax2.plot([], [], color='purple', ls="dashed", lw=1, label='$\sigma$ (Maxwell #3)') 
line_stress_ss, = ax2.plot([], [], color='black', ls="dashed", lw=1, label='$\sigma$ (spring)')
ax2.legend(loc='upper right')

# テキスト描画
eq_text = r'd$\sigma_j$/d$t$ = - $\sigma_j$/$\tau_j$'
ax2.text(0.1, 0.9, eq_text, transform=ax2.transAxes)
res_text = r'$E_i$ = {0:.2f} MPa, $E_\infty$ = {1:.2f} MPa'.format(insMod/10**6, infMod/10**6)
ax2.text(0.1, 0.8, res_text, transform=ax2.transAxes)

# テーブルの追加（グラフの下に配置）
param_table = ax2.table(
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

# アニメーション更新関数
def animate(i):
    # 歪みデータの更新
    line_strain.set_data(t[:i], e[:i])
    line_strain_e1.set_data(t[:i], e1[:i])
    line_strain_e2.set_data(t[:i], e2[:i])
    line_strain_e3.set_data(t[:i], e3[:i])

    # 応力データの更新
    line_stress.set_data(t[:i], s[:i])
    line_stress_s1.set_data(t[:i], s1[:i])
    line_stress_s2.set_data(t[:i], s2[:i])
    line_stress_s3.set_data(t[:i], s3[:i])
    line_stress_ss.set_data(t[:i], ss[:i])

    return line_strain, line_strain_e1, line_strain_e2, line_strain_e3, line_stress, line_stress_s1, line_stress_s2, line_stress_s3, line_stress_ss

# アニメーション実行   
ani = animation.FuncAnimation(fig, animate, frames=steps, interval=interval_ms, blit=True)

savefile = './mp4/genMaxwell_stepStrain.mp4'
ani.save(savefile, writer='ffmpeg', fps=fps, extra_args=['-r', '30'])

plt.tight_layout()
plt.show()