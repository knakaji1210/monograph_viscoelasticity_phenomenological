# Mittag=Leffler関数の表示

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import special
# import differint.differint as df # 古いバージョン
from differintP.functions import MittagLeffler

cmap = plt.get_cmap('Pastel1')

if __name__=='__main__':
    func_list = [r'$e^t$', r'$\cosh(t)$', r'$\sinh(t) / t$', r'$e^t (1 + \mathrm{erf}(\sqrt{t}))$']  # 関数のリスト
    alpha = [1.0, 2.0, 2.0, 0.5]  # Mittag-Leffler関数のパラメータ alpha のリスト
    beta = [1.0, 1.0, 2.0, 1.0]  # Mittag-Leffler関数のパラメータ beta のリスト

    # 計算条件の設定
    ts = 10**(-5)   # 下限（0にすると0除算が起こる関数があるので）
    te = 3          # 上限
    n = 1000        # 分割数
    t_vals = np.linspace(ts, te, n)

    x = np.zeros((4, n))  # X軸のデータ点を格納する配列
    y_orig = np.zeros((4, n))  # 元の関数の値を格納する配列
    y_ML = np.zeros((4, n))    # Mittag-Leffler関数の値を格納する配列

    # exp(x)
    x[0] = t_vals  # exp(x)のためにt=xとする
    y_orig[0] = np.exp(t_vals)
    y_ML[0] = MittagLeffler(alpha[0], beta[0], x[0], num_terms=50)   

    # cosh(x)
    x[1] = t_vals**2  # cosh(x)のためにx=t^2とする
    y_orig[1] = np.cosh(t_vals)
    y_ML[1] = MittagLeffler(alpha[1], beta[1], x[1], num_terms=50)

    # sinh(x)/x
    x[2] = t_vals**2  # sinh(x)/xのためにx=t^2とする
    y_orig[2] = np.sinh(t_vals) / t_vals
    y_ML[2] = MittagLeffler(alpha[2], beta[2], x[2], num_terms=50)

    # error関数
    x[3] = np.sqrt(t_vals)    # error関数のためにx
    y_ML[3] = MittagLeffler(alpha[3], beta[3], x[3], num_terms=50)
    y_orig[3] = np.exp(t_vals) * (1 + special.erf(np.sqrt(t_vals)))
        
    # グラフの描画
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.grid()
    ax.set_xlim(0, te)
    ax.set_ylim(0, 10)
    for i in range(4):
        ax.plot(t_vals, y_ML[i], label=f"Mittag-Leffler Function, $\\alpha={alpha[i]}$, $\\beta={beta[i]}$", color="gray", linestyle="-")
        ax.plot(t_vals, y_orig[i], label=f"Equivalent Function, {func_list[i]}", color=cmap(i/4), linestyle="--")
    ax.set_xlabel(r"$t$", fontsize=12)
    ax.set_ylabel(r"$E_{\alpha, \beta}(t)$", fontsize=12)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(fontsize=11)

    savefile = './png/MittagLeffler_multi.png'

    fig.savefig(savefile, dpi=300)

    # 表示
    plt.show()
