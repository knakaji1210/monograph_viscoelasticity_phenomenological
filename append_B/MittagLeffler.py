# Mittag=Leffler関数の表示

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import special
# import differint.differint as df # 古いバージョン
from differintP.functions import MittagLeffler

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
    select_text = 'Selection (exp: 0, cosh: 1, hypersinc: 2, error: 3, cos: 4, sin: 5): '
    try:
        select = int(input(select_text))
    except ValueError:
        select = 0

    func_list = [r'$e^t$', r'$\cosh(t)$', r'$\sinh(t) / t$', r'$e^t (1 + \mathrm{erf}(\sqrt{t}))$', r'$\cos(t)$', r'$\sin(t)$']  # 関数のリスト

    # 計算条件の設定
    ts = 10**(-5)   # 下限（0にすると0除算が起こる関数があるので）
    if select in [0, 1, 2, 3]:  # exp, cosh, sinh/t, errorの場合
        te = 3      # 上限
    else:
        te = 6.28   # 上限
    n = 1000        # 分割数
    t_vals = np.linspace(ts, te, n)
    
    if select == 0:     # exp(x)の場合
        alpha = 1.0  # Mittag-Leffler関数のパラメータ
        beta = 1.0   # Mittag-Leffler関数のパラメータ
        # X軸のデータ点を生成        
        x_vals = t_vals  # exp(x)のためにt=xとする        
        # 各計算の実行
        y_orig = np.exp(t_vals)
        y_ML = MittagLeffler(alpha, beta, x_vals, num_terms=50)
        savefile = './png/MittagLeffler_exp.png'
    elif select == 1:   # cosh(x)の場合
        alpha = 2.0  # Mittag-Leffler関数のパラメータ
        beta = 1.0   # Mittag-Leffler関数のパラメータ
        # X軸のデータ点を生成
        x_vals = t_vals**2  # cosh(x)のためにx=t^2とする
        # 各計算の実行
        y_orig = np.cosh(t_vals)
        y_ML = MittagLeffler(alpha, beta, x_vals, num_terms=50)
        savefile = './png/MittagLeffler_cosh.png'
    elif select == 2:   # sinh(x)/xの場合
        alpha = 2.0  # Mittag-Leffler関数のパラメータ
        beta = 2.0   # Mittag-Leffler関数のパラメータ
        # X軸のデータ点を生成
        x_vals = t_vals**2  # sinh(x)/xのためにx=t^2とする
        # 各計算の実行
        y_orig = np.sinh(t_vals) / t_vals
        y_ML = MittagLeffler(alpha, beta, x_vals, num_terms=50)
        savefile = './png/MittagLeffler_hypersinc.png'
    elif select == 3:   # error関数の場合
        alpha = 0.5  # Mittag-Leffler関数のパラメータ
        beta = 1.0   # Mittag-Leffler関数のパラメータ
        # X軸のデータ点を生成
        x_vals = np.sqrt(t_vals)    # error関数のためにx=sqrt(t)とする
        # 各計算の実行
        y_orig = np.exp(t_vals) * (1 + special.erf(np.sqrt(t_vals)))
        y_ML = MittagLeffler(alpha, beta, x_vals, num_terms=50)
        savefile = './png/MittagLeffler_error.png'
    elif select == 4:   # cos(x)の場合
        alpha = 2.0  # Mittag-Leffler関数のパラメータ
        beta = 1.0   # Mittag-Leffler関数のパラメータ
        # X軸のデータ点を生成
        x_vals = -t_vals**2  # cos(x)のためにx=-t^2とする
        # 各計算の実行
        # 組み込み関数で計算（ignore_special_cases=Trueを指定、指定しないとエラーが出る）
        y_orig = np.cos(t_vals)
        y_ML = MittagLeffler(alpha, beta, x_vals, num_terms=50, ignore_special_cases=True)
#        y_ML = MittagLeffler_self(x_vals, terms=50)  # 自作関数（MittagLeffler_self）で計算
        savefile = './png/MittagLeffler_cos.png'
    elif select == 5:   # sin(x)の場合
        alpha = 2.0  # Mittag-Leffler関数のパラメータ
        beta = 2.0   # Mittag-Leffler関数のパラメータ
        # X軸のデータ点を生成
        x_vals = -t_vals**2  # sin(x)のためにx=-t^2とする
        # 各計算の実行
        y_orig = np.sin(t_vals)
        y_ML = t_vals * MittagLeffler(alpha, beta, x_vals, num_terms=50, ignore_special_cases=True)
        savefile = './png/MittagLeffler_sin.png'
    else:
        pass

    ymax = np.max(y_orig)*1.1
    if select in [0, 1, 2, 3]:
        ylim = np.array([0, ymax])
    else:
        ylim = np.array([-ymax, ymax])

    # グラフの描画
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title(r'Mittag-Leffler Function, $E_{\alpha, \beta}(t)$', fontsize=14)
    ax.grid()
    ax.set_xlim(0, te)
    ax.set_ylim(ylim[0], ylim[1])
    ax.plot(t_vals, y_ML, label=f"Mittag-Leffler Function, $\\alpha={alpha}$, $\\beta={beta}$", color="gray", linestyle="-")
    ax.plot(t_vals, y_orig, label=f"Equivalent Function, $f(t)$ = {func_list[select]}", color="greenyellow", linestyle="--")
    ax.set_xlabel(r"$t$", fontsize=12)
    ax.set_ylabel(r"$E_{\alpha, \beta}(t)$", fontsize=12)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(fontsize=11)

    fig.savefig(savefile, dpi=300)

    # 表示
    plt.show()
