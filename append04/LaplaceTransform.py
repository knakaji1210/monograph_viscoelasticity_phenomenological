# Laplace Transform for various functions

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
 
# Laplace Transform
def calc_laplace(f_t, t, s):
# laplace_transformは (変換後の式, 収束域, 条件) を返すため、[0]で式のみを取得
    F_s = sp.laplace_transform(f_t, t, s)[0]
    return F_s

if __name__=='__main__':

    # 記号と関数の定義
#    t, s = sp.symbols('t s', real=True, positive=True)
# ヘビサイド関数をかける形に変更したため、次の段の表現に変更
    t, s = sp.symbols('t s', real=True)

    select_text = 'Selection (sin: 0, cos: 1, Heaviside: 2, ramp: 3, power: 4, exponential: 5, exponential*sine: 6, Dirac delta: 7): '  
    try:
        select = int(input(select_text))
    except ValueError:
        select = 0

    if select == 0:
        f_t = sp.sin(t) * sp.Heaviside(t)   # sine
        legend = 'sine'
        savefile = './png/Laplace_sine.png'
    if select == 1:
        f_t = sp.cos(t) * sp.Heaviside(t)   # cosine
        legend = 'cosine'
        savefile = './png/Laplace_cosine.png'
    if select == 2:
        f_t = sp.Heaviside(t)               # Heaviside
        legend = 'Heaviside'
        savefile = './png/Laplace_Heaviside.png'
    if select == 3:
        f_t = t * sp.Heaviside(t)           # ramp
        legend = 'ramp'
        savefile = './png/Laplace_ramp.png'
    if select == 4:
        f_t = t**2 * sp.Heaviside(t)        # power
        legend = 'power'
        savefile = './png/Laplace_power.png'
    if select == 5:
        f_t = sp.exp(-t) * sp.Heaviside(t)  # exponential
        legend = 'exponential'
        savefile = './png/Laplace_exponential.png'
    if select == 6:
        f_t = sp.exp(-t) * sp.sin(t) * sp.Heaviside(t)
        legend = 'exponential * sine'
        savefile = './png/Laplace_exp_sine.png'
    if select == 7:
        f_t = sp.DiracDelta(t)              # Dirac delta
        legend = 'Dirac delta'
        savefile = './png/Laplace_Dirac_delta.png'

    # ラプラス変換の実行
    F_s = calc_laplace(f_t, t, s)

    # グラフ描画用のデータ生成
    t_vals = np.linspace(-1, 10, 400)
    s_vals = np.linspace(10**(-3), 10, 400)  # s=0での発散を避けるため10^(-3)から開始

    # SymPyの数式をNumPyで計算できる関数に変換 (lambdify)
    if select == 7: # Dirac deltaの場合のみ例外的に処理する（numpyでDirac deltaを扱えないため）
        f_y = np.zeros_like(t_vals)
        F_y = np.full_like(s_vals, float(F_s))
    else:
        f_num = sp.lambdify(t, f_t, 'numpy')
        f_y = f_num(t_vals)
        F_num = sp.lambdify(s, F_s, 'numpy')
        F_y = F_num(s_vals)

    xlim = 10
    # Dirac deltaの場合はy軸の上限を1.2に設定
    if select == 7:
        y1lim = 1.2
    else:
        y1lim = np.max(np.abs(f_y)) * 1.2
    # Heatvisideやramp、powerはy軸の値が大きくなるため、上限を10に設定
    if select == 2 or select == 3 or select == 4:
        y2lim = 10
    else:
        y2lim = np.max(np.abs(F_y)) * 1.2

    fig = plt.figure(figsize=(12,4), tight_layout=True)

    ax1 = fig.add_subplot(121)
#    ax1.set_title(r'Original function: $f(t)$')
    ax1.set_xlabel(r'$t$ /s')
    ax1.set_ylabel(r'$f(t)$')
    ax1.set_xlim(-1, xlim)
    ax1.set_ylim(-y1lim, y1lim)
    ax1.plot(t_vals, f_y, label=legend, color='blue')
    ax1.legend(loc='upper right')
    ax1.hlines(0, -1, xlim, color='black', lw=0.5, ls='dashed')
    ax1.vlines(0, -y1lim, y1lim, color='black', lw=0.5, ls='dashed')
    if select == 7:  # Dirac deltaの場合はt=0に矢印を描画
        ax1.vlines(0, 0, y1lim, color='blue')

    ax2 = fig.add_subplot(122)
#    ax2.set_title(r'Laplace transform: $F(s)$')
    ax2.set_xlabel(r'$s$')
    ax2.set_ylabel(r'$F(s)$')
    ax2.set_xlim(-1, xlim)
    ax2.set_ylim(-y2lim, y2lim)
    ax2.plot(s_vals, F_y, label='Laplace', color='red')
    ax2.legend(loc='upper right')
    ax2.hlines(0, -1, xlim, color='black', lw=0.5, ls='dashed')
    ax2.vlines(0, -y2lim, y2lim, color='black', lw=0.5, ls='dashed')
 
    fig.savefig(savefile, dpi=300)

    plt.show()