# RL微分とCaputo微分の比較

'''
differintPでcaputo微分はポイントでしか行えない
そこでRL微分もポイントで行い、Caputo微分と比較する
'''

import numpy as np
import matplotlib.pyplot as plt
from differintP import RLpoint        # Riemann-Liouville微分
from differintP import CaputoL1point  # Caputo微分（L1）

# 対象とする関数
def f0(t):
    return t**0

def f1(t):
    return t

def f2(t):
    return t**2

def f3(t):
    return np.sqrt(t)

def f4(t):
    return np.exp(t)

def f5(t):
    return np.sin(t)

def f6(t):
    return np.cos(t)

def f7(t):
    return t**2 + 0.5

# RL微分の計算
def calc_RL_derivative(func, nu, t_vals):
    RL_derivs = []
    for t in t_vals:
        if t == t_vals[0]:
            val = 0.0
        else:
            val = RLpoint(nu, func, t_vals[0], t, len(t_vals))
        RL_derivs.append(val)
    return RL_derivs

# Caputo微分の計算
def calc_Caputo_derivative(func, nu, t_vals):
    Caputo_derivs = []
    for t in t_vals:
        if t == t_vals[0]:
            val = 0.0
        else:
            val = CaputoL1point(nu, func, t_vals[0], t, len(t_vals))
        Caputo_derivs.append(val)
    return Caputo_derivs

if __name__=='__main__':
    select_text = 'Selection (1: 0, t: 1, t^2: 2, sqrt: 3, exp: 4, sin: 5, cos: 6, t^2 + 0.5: 7): '  
    try:
        select = int(input(select_text))
    except ValueError:
        select = 0

    try:
        nu = float(input("Enter the order of the derivative (0 < nu < 1): "))
        if not (0 < nu < 1):
            raise ValueError("Order must be between 0 and 1.")
    except ValueError:
        nu = 0.5  # Default value if input is invalid

    # 時間領域と関数の定義
    ts = 0.0
    if select in [0, 1, 2, 3]:
        te = 3.0  # 上限
    elif select in [5, 6]:  # sin, cosの場合
        te = 6.28  # 上限
    else:
        te = 1.0  # 上限
    n = 1000
    t_vals = np.linspace(ts, te, n)

    func_list = [r'$1$', r'$t$', r'$t^2$', r'$\sqrt{t}$', r'$e^t$', r'$\sin(t)$', r'$\cos(t)$', r'$t^2$ + 0.5']  # 関数のリスト

    # 各計算の実行
    if select == 0:
        y_orig = f0(t_vals)
        y_RL = calc_RL_derivative(f0, nu, t_vals)
        y_Caputo = calc_Caputo_derivative(f0, nu, t_vals)
        ylim = np.array([-0.1, 3.0])
        savefile = './png/RL_Caputo_const.png'
    elif select == 1:
        y_orig = f1(t_vals)
        y_RL = calc_RL_derivative(f1, nu, t_vals)
        y_Caputo = calc_Caputo_derivative(f1, nu, t_vals)
        ylim = np.array([-0.1, 3.0])
        savefile = './png/RL_Caputo_linear.png'
    elif select == 2:
        y_orig = f2(t_vals)
        y_RL = calc_RL_derivative(f2, nu, t_vals)
        y_Caputo = calc_Caputo_derivative(f2, nu, t_vals)
        ylim = np.array([-0.1, 10.0])
        savefile = './png/RL_Caputo_square.png'
    elif select == 3:
        y_orig = f3(t_vals)
        y_RL = calc_RL_derivative(f3, nu, t_vals)
        y_Caputo = calc_Caputo_derivative(f3, nu, t_vals)
        ylim = np.array([-0.1, 3.0])
        savefile = './png/RL_Caputo_sqrt.png'
    elif select == 4:
        y_orig = f4(t_vals)
        y_RL = calc_RL_derivative(f4, nu, t_vals)
        y_Caputo = calc_Caputo_derivative(f4, nu, t_vals)
        ylim = np.array([-0.1, 10.0])
        savefile = './png/RL_Caputo_exp.png'
    elif select == 5:
        y_orig = f5(t_vals)
        y_RL = calc_RL_derivative(f5, nu, t_vals)
        y_Caputo = calc_Caputo_derivative(f5, nu, t_vals)
        ylim = np.array([-2.0, 2.0])
        savefile = './png/RL_Caputo_sin.png'
    elif select == 6:
        y_orig = f6(t_vals)
        y_RL = calc_RL_derivative(f6, nu, t_vals)
        y_Caputo = calc_Caputo_derivative(f6, nu, t_vals)
        ylim = np.array([-2.0, 2.0])
        savefile = './png/RL_Caputo_cos.png'
    elif select == 7:
        y_orig = f7(t_vals)
        y_RL = calc_RL_derivative(f7, nu, t_vals)
        y_Caputo = calc_Caputo_derivative(f7, nu, t_vals)
        ylim = np.array([-0.1, 3.0])
        savefile = './png/RL_Caputo_square_const.png'

    # グラフの描画
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.grid()

    ax.plot(t_vals, y_orig, label="Original: $f(t)$ = {}".format(func_list[select]), color="black", linestyle="--")
    ax.plot(t_vals, y_RL, label=f"{nu}-order Derivative (RL)", color="greenyellow", linewidth=2)
    ax.plot(t_vals, y_Caputo, label=f"{nu}-order Derivative (Caputo)", color="orange", linewidth=2, linestyle="-.")

    ax.set_title(f"Fractional Derivatives of $f(t)$", fontsize=14)
    ax.set_xlabel(r"$t$", fontsize=12)
    ax.set_ylabel(r"$D^{\nu} f(t)$", fontsize=12)
    ax.set_ylim(ylim[0], ylim[1])
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(fontsize=11)

    fig.savefig(savefile, dpi=300)

    # 表示
    plt.show()
