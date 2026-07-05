# Mittag=Leffler関数の漸近的振る舞い

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import special
# import differint.differint as df # 古いバージョン
from differintP.functions import MittagLeffler

if __name__=='__main__':
    select_text = 'Selection (linear: 0, log: 1): '  
    try:
        select = int(input(select_text))
    except ValueError:
        select = 0

    try:
        a = float(input('Input the scale factor (alpha = 1.0): '))
    except ValueError:
        a = 1.0  # Default value if input is invalid

    try:
        nu = float(input('Input the order of fractional derivative (0 < nu < 1): '))
    except ValueError:
        nu = 0.5  # Default value if input is invalid

    func_list = [r'exp($-{\alpha t}$)', r'$E_{\nu, 1}(-\alpha t^{\nu})$', r'exp($-{\alpha t^{\nu}}$)', r'$-\alpha t^{\nu}$']  # 関数のリスト

    if select == 0:  # 線形スケールの場合
        ts = 0      # 下限
        te = 3      # 上限
        n = 1000    # 分割数
        t_vals = np.linspace(ts, te, n)
    else:  # 対数スケールの場合
        log10_ts = -1.2   # 下限（0にすると0除算が起こる関数があるので）
        log10_te = 1.2    # 上限
        n = 1000        # 分割数
        t_vals = np.logspace(log10_ts, log10_te, n)
        log10_t_vals = np.log10(t_vals)

    '''
    ここで用いてる differintP.functions.MittagLeffler()は負の引数に対して
    tの値が大きくなると、計算が不安定になるので、範囲を狭めて表示している。
    ./legacy/MittagLeffler_asymptotic_GoogleAI.pyでは小さい負の値に対しては級数展開、
    大きい負の値に対しては漸近展開を用いて計算しているのでより広い範囲で表示できる。
    '''

    x_vals = -a * t_vals**nu

    y_exp = np.exp(-a*t_vals)
    log10_y_exp = np.log10(y_exp)
    y_ML = MittagLeffler(nu, 1, x_vals, num_terms=100)
    log10_y_ML = np.log10(y_ML)
    y_KWW = np.exp(-a*t_vals**nu)
    log10_y_KWW = np.log10(y_KWW)
    y_power = 1.0 / (a * special.gamma(1 - nu)) * t_vals**(-nu)
    log10_y_power = np.log10(y_power)


    # グラフの描画
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title(f'Asymptotic Behavior of Mittag-Leffler Function ($\\alpha={a}, \\nu={nu}$)', fontsize=14)
    ax.grid()
    if select == 0:
        ax.set_xlim(ts, te)
        ax.set_ylim(0, 1.2)
        ax.plot(t_vals, y_exp, label=f"{func_list[0]}", color="gray", linestyle="-")
        ax.plot(t_vals, y_ML, label=f"{func_list[1]}", color="greenyellow", linestyle="-")
        ax.plot(t_vals, y_KWW, label=f"{func_list[2]}", color="orange", linestyle="--")
        ax.plot(t_vals, y_power, label=f"{func_list[3]}", color="cyan", linestyle="--")
        ax.set_xlabel(r"$t$", fontsize=12)
        ax.set_ylabel(r"$f(t)$", fontsize=12)
        ax.grid(True, linestyle=":", alpha=0.6)
        ax.legend(fontsize=11)
        savefile = './png/MittagLeffler_asymptotic_linear.png'

    else:
        ax.set_xlim(log10_ts, log10_te)
        ax.set_ylim(-1, 0.5)
        ax.plot(log10_t_vals, log10_y_exp, label=f"{func_list[0]}", color="gray", linestyle="-")
        ax.plot(log10_t_vals, log10_y_ML, label=f"{func_list[1]}", color="greenyellow", linestyle="-")
        ax.plot(log10_t_vals, log10_y_KWW, label=f"{func_list[2]}", color="orange", linestyle="--")
        ax.plot(log10_t_vals, log10_y_power, label=f"{func_list[3]}", color="cyan", linestyle="--")
        ax.set_xlabel(r"log$_{10} (t)$", fontsize=12)
        ax.set_ylabel(r"log$_{10} (f(t))$", fontsize=12)
        savefile = './png/MittagLeffler_asymptotic_log.png'

    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(fontsize=11)

    fig.savefig(savefile, dpi=300)

    # 表示
    plt.show()
