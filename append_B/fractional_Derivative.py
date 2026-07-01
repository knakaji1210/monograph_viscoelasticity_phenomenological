# 非整数階微分の計算とグラフ化（Riemann-Liouville定義）

'''
differintP：非整数階微積分（fractional calculus）を計算するためのPythonライブラリ
（differintは古いバージョンでアップデートがかからないらしいのでこちらに移行した）
alpha > 0 の場合は非整数階微分、alpha < 0 の場合は非整数階積分を計算

df.GL / df.GLI：Grünwald-Letnikov定義、alphaは全ての実数値に対応
df.RL / df.RLI：Riemann-Liouville定義、alpha < 1
df.CaputoL1point：Caputo（L1）定義、0 < alpha < 1
df.CaputoL2point：Caputo（L2）定義、1 < alpha < 2
df.CaputoL2Cpoint：Caputo（L2C）定義、0 < alpha < 1, 1 < alpha < 2

ここではdf.RLを利用
なお、Caputoの実装は不完全のようだ
'''

import numpy as np
import matplotlib.pyplot as plt
import differintP as df

# 対象の関数
def f0(t):
    return t**0
# differintの場合は下でも大丈夫だったが、differintPではエラーになるので、上のように変更
#    return 1

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

cmap = plt.get_cmap('coolwarm')

if __name__=='__main__':
    select_text = 'Selection (1: 0, t: 1, t^2: 2, sqrt: 3, exp: 4, sin: 5, cos: 6): '  
    try:
        select = int(input(select_text))
    except ValueError:
        select = 0

    # 階数の配列
    nu_arr = np.array([0.2, 0.4, 0.6, 0.8, 1.0])
    # 古いバージョンでは微分の方は1.0だ計算できず、代わりに0.99で計算していたが、differintPでは1.0でも計算できるようになった
    func_list = [r'$1$', r'$t$', r'$t^2$', r'$\sqrt{t}$', r'$e^t$', r'$\sin(t)$', r'$\cos(t)$']  # 関数のリスト

    # X軸のデータ点を生成
    ts = 0         # 下限
    if select == 4:  # exp(t)の場合は上限を小さくする
        te = 2.0      # 上限
    elif select == 5 or select == 6:  # sin(t), cos(t)の場合
        te = 6.28      # 上限
    else:    
        te = 4.0      # 上限
    n = 1000        # 分割数
    t_vals = np.linspace(ts, te, n)

    # 各計算の実行
    if select == 0:
        # 元の関数の値
        y_orig = np.ones_like(t_vals)
        y_int_array = np.zeros((len(nu_arr), n))  # 階数ごとの積分結果を格納する配列
        for i in range(len(nu_arr)):
            y_int_array[i] = df.RL(nu_arr[i], f0, ts, te, n)
        ylim = np.array([-0.1, 3.0])
        savefile = './png/fractional_derivative_const.png'
    elif select == 1:
        y_orig = f1(t_vals)
        y_int_array = np.zeros((len(nu_arr), n))  # 階数ごとの積分結果を格納する配列
        for i in range(len(nu_arr)):
            y_int_array[i] = df.RL(nu_arr[i], f1, ts, te, n)
        ylim = np.array([-0.1, 3.0])
        savefile = './png/fractional_derivative_linear.png'
    elif select == 2:
        y_orig = f2(t_vals)
        y_int_array = np.zeros((len(nu_arr), n))  # 階数ごとの積分結果を格納する配列
        for i in range(len(nu_arr)):
            y_int_array[i] = df.RL(nu_arr[i], f2, ts, te, n)
        ylim = np.array([-1.0, 10.0])
        savefile = './png/fractional_derivative_square.png'
    elif select == 3:
        y_orig = f3(t_vals)
        y_int_array = np.zeros((len(nu_arr), n))  # 階数ごとの積分結果を格納する配列
        for i in range(len(nu_arr)):
            y_int_array[i] = df.RL(nu_arr[i], f3, ts, te, n)
        ylim = np.array([-0.1, 3.0])
        savefile = './png/fractional_derivative_sqrt.png'
    elif select == 4:
        y_orig = f4(t_vals)
        y_int_array = np.zeros((len(nu_arr), n))  # 階数ごとの積分結果を格納する配列
        for i in range(len(nu_arr)):
            y_int_array[i] = df.RL(nu_arr[i], f4, ts, te, n)
        ylim = np.array([-0.1, 10.0])
        savefile = './png/fractional_derivative_exp.png'
    elif select == 5:
        y_orig = f5(t_vals)
        y_int_array = np.zeros((len(nu_arr), n))  # 階数ごとの積分結果を格納する配列
        for i in range(len(nu_arr)):
            y_int_array[i] = df.RL(nu_arr[i], f5, ts, te, n)
        ylim = np.array([-2.0, 2.0])
        savefile = './png/fractional_derivative_sin.png'
    elif select == 6:
        y_orig = f6(t_vals)
        y_int_array = np.zeros((len(nu_arr), n))  # 階数ごとの積分結果を格納する配列
        for i in range(len(nu_arr)):
            y_int_array[i] = df.RL(nu_arr[i], f6, ts, te, n)
        ylim = np.array([-2.0, 2.0])
        savefile = './png/fractional_derivative_cos.png'
    else:
        pass

    # グラフの描画
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.grid()

    ax.plot(t_vals, y_orig, label="Original: $f(t)$ = {}".format(func_list[select]), color="black", linestyle="--", zorder=2)
    for i in range(len(nu_arr)):
        ax.plot(t_vals, y_int_array[i], label=f"{nu_arr[i]}-order Derivative (RL)", color=cmap(i/len(nu_arr)), linewidth=2, zorder=1)
    ax.set_title(f"Fractional Derivative of $f(t)$", fontsize=14)
    ax.set_xlabel(r"$t$", fontsize=12)
    ax.set_ylabel(r"$D^{\nu} f(t)$", fontsize=12)
    ax.set_ylim(ylim[0], ylim[1])
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(fontsize=11)

    fig.savefig(savefile, dpi=300)

    # 表示
    plt.show()
