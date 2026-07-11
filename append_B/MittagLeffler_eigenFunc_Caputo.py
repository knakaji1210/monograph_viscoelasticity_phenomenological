# 非整数階微分の固有関数としてのMittag=Leffler関数の表示（Caputo微分（L1））

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import special
# import differint.differint as df # 古いバージョン
from differintP.functions import MittagLeffler

cmap = plt.get_cmap('coolwarm')

if __name__=='__main__':
    select_text = 'Selection (linear: 0, log: 1): '  
    try:
        select = int(input(select_text))
    except ValueError:
        select = 0
 
    try:
        a = float(input('Input the scale factor a (a > 0): '))
    except ValueError:
        a = 1.0

    nu_arr = np.array([0.05, 0.2, 0.4, 0.6, 0.8, 0.99])  # 階数の配列 
    func_list = [r'$e^{\alpha t}$', r'$E_{\nu, 1}(\alpha t^{\nu})$']  # 関数のリスト

    if select == 0:  # 線形スケールの場合
        ts = 0      # 下限
        te = 3      # 上限
        n = 1000    # 分割数
        t_vals = np.linspace(ts, te, n)
    elif select == 1:  # 対数スケールの場合
        log10_ts = -3   # 下限（0にすると0除算が起こる関数があるので）
        log10_te = 3    # 上限
        n = 1000        # 分割数
        t_vals = np.logspace(log10_ts, log10_te, n)
        log10_t_vals = np.log10(t_vals)

    y_exp = np.exp(a * t_vals)
    log10_y_exp = np.log10(y_exp)
    y_eiginFunc = np.zeros((len(nu_arr), len(t_vals)))  # 階数ごとの積分結果を格納する配列
    for i in range(len(nu_arr)):
        nu = nu_arr[i]
        x_vals = a * t_vals**nu
        y_eiginFunc[i] = MittagLeffler(nu, 1, x_vals, num_terms=50)
    log10_y_eiginFunc = np.log10(y_eiginFunc)
    ymax_linear = np.max(y_eiginFunc[-1])*1.1
    ymin_log = np.min(log10_y_eiginFunc[0])-3.0
    ymax_log = np.max(log10_y_eiginFunc[0])*2.0

# グラフの描画
fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111)
ax.set_title(f'Eigen Function of $D^{{1}}$ and $D^{{\\nu}}_{{C}}$ (Caputo)', fontsize=14)
ax.grid()
if select == 0:
    ax.set_xlim(ts, te)
    ax.set_ylim(0, ymax_linear)
    ax.plot(t_vals, y_exp, label=f"eigen Function of $D^{{1}} f(t)$, {func_list[0]}", color="gray", linestyle="-")
    for i in range(len(nu_arr)):
        nu = nu_arr[i]
        ax.plot(t_vals, y_eiginFunc[i], label=f"eigen Function of $D^{{\\nu}}_{{C}} f(t), \\nu={nu}$", color=cmap(i/len(nu_arr)), linestyle="--")
    ax.set_xlabel(r"$t$", fontsize=12)
    ax.set_ylabel(f"{func_list[0]}, {func_list[1]}", fontsize=12)
    savefile = './png/MittagLeffler_eigenFunc_Caputo_linear.png'
elif select == 1:
    ax.set_xlim(log10_ts, log10_te)
    ax.set_ylim(ymin_log, ymax_log)
    ax.plot(log10_t_vals, log10_y_exp, label=f"eigen Function of $D^{{1}} f(t)$, {func_list[0]}", color="gray", linestyle="-")
    for i in range(len(nu_arr)):
        nu = nu_arr[i]
        ax.plot(log10_t_vals, log10_y_eiginFunc[i], label=f"eigen Function of $D^{{\\nu}}_{{C}} f(t), \\nu={nu}$", color=cmap(i/len(nu_arr)), linestyle="--")
    ax.set_xlabel(r"log$_{10} (t)$", fontsize=12)
    ax.set_ylabel(r"log$_{10} ($" + f"{func_list[0]}, {func_list[1]})", fontsize=12)
#    ax.set_ylabel(rf"log$_{{10}}$ ({func_list[0]}, {func_list[1]})$", fontsize=12)
    savefile = './png/MittagLeffler_eigenFunc_Caputo_log.png'
ax.grid(True, linestyle=":", alpha=0.6)
ax.legend(fontsize=11)

fig.savefig(savefile, dpi=300)

# 表示
plt.show()
