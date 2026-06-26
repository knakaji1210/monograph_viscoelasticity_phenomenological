# Mittag=Leffler関数の表示

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import special
import differint.differint as df


# 計算条件の設定
a = 10**(-3)  # 下限（0にすると0除算が起こる関数があるので）
b = 3         # 上限
n = 500       # 分割数

if __name__=='__main__':
    # X軸のデータ点を生成
    t_vals = np.linspace(a, b, n)
    # exp(x)
    x0_vals = t_vals  # exp(x)のためにt=xとする        
    y0 = df.MittagLeffler(1.0, 1.0, x0_vals, num_terms=50)
    y0_orig = np.exp(t_vals)
    f0_label = r"$f(t) = e^t$"

    # cosh(x)
    x1_vals = t_vals**2  # cosh(x)のためにx=t^2とする
    y1 = df.MittagLeffler(2.0, 1.0, x1_vals, num_terms=50)
    y1_orig = np.cosh(t_vals)
    f1_label = r"$f(t) = \cosh(t)$"

    # sinh(x)/x
    x2_vals = t_vals**2  # sinh(x)/xのためにx=t^2とする
    y2 = df.MittagLeffler(2.0, 2.0, x2_vals, num_terms=50)
    y2_orig = np.sinh(t_vals) / t_vals
    f2_label = r"$f(t) = \sinh(t) / t$"

    # error関数
    x3_vals = np.sqrt(t_vals)    # error関数のためにx
    y3 = df.MittagLeffler(0.5, 1.0, x3_vals, num_terms=50)
    y3_orig = np.exp(t_vals) * (1 + special.erf(np.sqrt(t_vals)))
    f3_label = r"$f(t) = e^t (1 + \mathrm{erf}(\sqrt{t}))$"
        

    # グラフの描画
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.grid()
    ax.set_xlim(0, b)
    ax.set_ylim(0, 10)
    ax.plot(t_vals, y0, label=f"Mittag-Leffler Function, $\\alpha=1.0$, $\\beta=1.0$", color="gray", linestyle="-")
    ax.plot(t_vals, y1, label=f"Mittag-Leffler Function, $\\alpha=2.0$, $\\beta=1.0$", color="gray", linestyle="-")
    ax.plot(t_vals, y2, label=f"Mittag-Leffler Function, $\\alpha=2.0$, $\\beta=2.0$", color="gray", linestyle="-")
    ax.plot(t_vals, y3, label=f"Mittag-Leffler Function, $\\alpha=0.5$, $\\beta=1.0$", color="gray", linestyle="-")
    ax.plot(t_vals, y0_orig, label=f"Equivalent Function, {f0_label}", color="pink", linestyle="--")
    ax.plot(t_vals, y1_orig, label=f"Equivalent Function, {f1_label}", color="aqua", linestyle="--")
    ax.plot(t_vals, y2_orig, label=f"Equivalent Function, {f2_label}", color="greenyellow", linestyle="--")
    ax.plot(t_vals, y3_orig, label=f"Equivalent Function, {f3_label}", color="orange", linestyle="--")
    ax.set_xlabel(r"$t$", fontsize=12)
    ax.set_ylabel(r"$E_{\alpha, \beta}(t)$", fontsize=12)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(fontsize=11)

    savefile = './png/MittagLeffler_multi.png'

    fig.savefig(savefile, dpi=300)

    # 表示
    plt.show()
