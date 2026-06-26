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

# Mittag-Leffler関数の定義（cos(x)が組み込みの関数でエラーが出たので、試しに自作）
# 組み込み関数でも「ignore_special_cases=True」を指定すればできることがわかった
# が、せっかくなので残す
def MittagLeffler_self(z, terms):
    """パラメータ alpha=2, beta=1 の Mittag-Leffler 関数を級数和で計算"""
    result = np.zeros_like(z, dtype=float)
    for k in range(terms):
        # 1 / Gamma(2k + 1) = 1 / (2k)!
        factorial = float(math.factorial(2 * k))
        result += (z**k) / factorial
    return result

if __name__=='__main__':
    select_text = 'Selection (exp(x): 0, cosh(x): 1, sinh(x)/x: 2, error: 3, cos(x): 4): '
    try:
        select = int(input(select_text))
    except ValueError:
        select = 0
    
    if select == 0:     # exp(x)の場合
        alpha = 1.0  # Mittag-Leffler関数のパラメータ
        beta = 1.0   # Mittag-Leffler関数のパラメータ
        # X軸のデータ点を生成
        t_vals = np.linspace(a, b, n)
        x_vals = t_vals  # exp(x)のためにt=xとする        
        # 各計算の実行
        y_ML = df.MittagLeffler(alpha, beta, x_vals, num_terms=50)
        y_orig = np.exp(t_vals)
        ymax = np.max(y_orig)*1.1
        ylim = np.array([0, ymax])
        func_label = r"$f(t) = e^t$"
        savefile = './png/MittagLeffler_exp.png'
    elif select == 1:   # cosh(x)の場合
        alpha = 2.0  # Mittag-Leffler関数のパラメータ
        beta = 1.0   # Mittag-Leffler関数のパラメータ
        # X軸のデータ点を生成
        t_vals = np.linspace(a, b, n)
        x_vals = t_vals**2  # cosh(x)のためにx=t^2とする
        # 各計算の実行
        y_ML = df.MittagLeffler(alpha, beta, x_vals, num_terms=50)
        y_orig = np.cosh(t_vals)
        ymax = np.max(y_orig)*1.1
        ylim = np.array([0, ymax])
        func_label = r"$f(t) = \cosh(t)$"
        savefile = './png/MittagLeffler_cosh.png'
    elif select == 2:   # sinh(x)/xの場合
        alpha = 2.0  # Mittag-Leffler関数のパラメータ
        beta = 2.0   # Mittag-Leffler関数のパラメータ
        # X軸のデータ点を生成
        t_vals = np.linspace(a, b, n)
        x_vals = t_vals**2  # sinh(x)/xのためにx=t^2とする
        # 各計算の実行
        y_ML = df.MittagLeffler(alpha, beta, x_vals, num_terms=50)
        y_orig = np.sinh(t_vals) / t_vals
        ymax = np.max(y_orig)*1.1
        ylim = np.array([0, ymax])
        func_label = r"$f(t) = \sinh(t) / t$"
        savefile = './png/MittagLeffler_hypersinc.png'
    elif select == 3:   # error関数の場合
        alpha = 0.5  # Mittag-Leffler関数のパラメータ
        beta = 1.0   # Mittag-Leffler関数のパラメータ
        # X軸のデータ点を生成
        t_vals = np.linspace(a, b, n)
        x_vals = np.sqrt(t_vals)    # error関数のためにx=sqrt(t)とする
        # 各計算の実行
        y_ML = df.MittagLeffler(alpha, beta, x_vals, num_terms=50)
        y_orig = np.exp(t_vals) * (1 + special.erf(np.sqrt(t_vals)))
        ymax = np.max(y_orig)*1.1
        ylim = np.array([0, ymax])
        func_label = r"$f(t) = e^t (1 + \mathrm{erf}(\sqrt{t}))$"
        savefile = './png/MittagLeffler_error.png'
    elif select == 4:   # cos(x)の場合
        alpha = 2.0  # Mittag-Leffler関数のパラメータ
        beta = 1.0   # Mittag-Leffler関数のパラメータ
        # X軸のデータ点を生成
        t_vals = np.linspace(a, b, n)
        x_vals = -t_vals**2  # cos(x)のためにx=-t^2とする
        # 各計算の実行
#        y_ML = df.MittagLeffler(alpha, beta, x_vals, num_terms=50)     # 組み込みだとエラーが出た
        y_ML = df.MittagLeffler(alpha, beta, x_vals, num_terms=50, ignore_special_cases=True)  # 組み込み関数で計算（ignore_special_cases=Trueを指定）できた
#        y_ML = MittagLeffler_self(x_vals, terms=50)  # 自作関数で計算
        y_orig = np.cos(t_vals)
        ymax = np.max(y_orig)*1.1
        ylim = np.array([-ymax, ymax])
        func_label = r"$f(t) = \cos(t)$"
        savefile = './png/MittagLeffler_cos.png'

    # グラフの描画
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.grid()
    ax.set_xlim(0, b)
    ax.set_ylim(ylim[0], ylim[1])
    ax.plot(t_vals, y_ML, label=f"Mittag-Leffler Function, $\\alpha={alpha}$, $\\beta={beta}$", color="gray", linestyle="-")
    ax.plot(t_vals, y_orig, label=f"Equivalent Function, {func_label}", color="greenyellow", linestyle="--")
    ax.set_xlabel(r"$t$", fontsize=12)
    ax.set_ylabel(r"$E_{\alpha, \beta}(t)$", fontsize=12)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(fontsize=11)

    fig.savefig(savefile, dpi=300)

    # 表示
    plt.show()
