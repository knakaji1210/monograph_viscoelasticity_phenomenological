# 非整数階積分の計算とグラフ化

'''
differint：非整数階微積分（fractional calculus）を計算するためのPythonライブラリ
alpha > 0 の場合は非整数階微分、alpha < 0 の場合は非整数階積分を計算

df.GL / df.GLI：Grünwald-Letnikov定義、alphaの制限なし
df.RL：Riemann-Liouville定義、、alphaの制限なし？
df.CaputoL1point：Caputo（L1）定義、0 < alpha < 1
df.CaputoL2point：Caputo（L2）定義、0 < alpha < 1
df.CaputoL2Cpoint：Caputo（L2C）定義、0 < alpha < 1, 1 < alpha < 2

ここではdf.RLを利用
'''

import numpy as np
import matplotlib.pyplot as plt
import differint.differint as df

# 対象の関数
def f0(t):
    return 1

def f1(t):
    return t

def f2(t):
    return t**2

def f3(t):
    return np.sqrt(t)

def f4(t):
    return np.exp(t)

def f5(t):
    return np.sin(2*np.pi*t/4)

# 計算条件の設定
a = 0         # 下限
b = 4         # 上限
n = 500       # 分割数

cmap = plt.get_cmap('coolwarm')

if __name__=='__main__':
    select_text = 'Selection (1: 0, x: 1, x^2: 2, sqrt(x): 3, exp(x): 4, sin(x): 5): '  
    try:
        select = int(input(select_text))
    except ValueError:
        select = 0

    lambda_arr = np.array([0.2, 0.4, 0.6, 0.8, 1.0])  # 階数の配列 
    func_list = [r'$1$', r'$x$', r'$x^2$', r'$\sqrt{x}$', r'$e^x$', r'$\sin(x)$']  # 関数のリスト

    # X軸のデータ点を生成
    t_vals = np.linspace(a, b, n)

    # 各計算の実行
    if select == 0:
        # 元の関数の値
        y_orig = np.ones_like(t_vals)
        y_int_array = np.zeros((len(lambda_arr), n))  # 階数ごとの積分結果を格納する配列
        for i in range(len(lambda_arr)):
            y_int_array[i] = df.RL(-lambda_arr[i], f0, a, b, n)
        savefile = './png/fractional_integral_const.png'
    elif select == 1:
        y_orig = f1(t_vals)
        y_int_array = np.zeros((len(lambda_arr), n))  # 階数ごとの積分結果を格納する配列
        for i in range(len(lambda_arr)):
            y_int_array[i] = df.RL(-lambda_arr[i], f1, a, b, n)
        savefile = './png/fractional_integral_linear.png'
    elif select == 2:
        y_orig = f2(t_vals)
        y_int_array = np.zeros((len(lambda_arr), n))  # 階数ごとの積分結果を格納する配列
        for i in range(len(lambda_arr)):
            y_int_array[i] = df.RL(-lambda_arr[i], f2, a, b, n)
        savefile = './png/fractional_integral_square.png'
    elif select == 3:
        y_orig = f3(t_vals)
        y_int_array = np.zeros((len(lambda_arr), n))  # 階数ごとの積分結果を格納する配列
        for i in range(len(lambda_arr)):
            y_int_array[i] = df.RL(-lambda_arr[i], f3, a, b, n)
        savefile = './png/fractional_integral_sqrt.png'
    elif select == 4:
        y_orig = f4(t_vals)
        y_int_array = np.zeros((len(lambda_arr), n))  # 階数ごとの積分結果を格納する配列
        for i in range(len(lambda_arr)):
            y_int_array[i] = df.RL(-lambda_arr[i], f4, a, b, n)
        savefile = './png/fractional_integral_exp.png'
    elif select == 5:
        y_orig = f5(t_vals)
        y_int_array = np.zeros((len(lambda_arr), n))  # 階数ごとの積分結果を格納する配列
        for i in range(len(lambda_arr)):
            y_int_array[i] = df.RL(-lambda_arr[i], f5, a, b, n)
        savefile = './png/fractional_integral_sin.png'
    else:
        pass

    # グラフの描画
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.grid()

    ax.plot(t_vals, y_orig, label="Original: $f(t)$ = {}".format(func_list[select]), color="black", linestyle="--")
    for i in range(len(lambda_arr)):
        ax.plot(t_vals, y_int_array[i], label=f"{lambda_arr[i]}-order Integral (RL)", color=cmap(i/len(lambda_arr)), linewidth=2)

    ax.set_title(f"Fractional Integral of $f(t)$", fontsize=14)
    ax.set_xlabel(r"$t$", fontsize=12)
    ax.set_ylabel(r"$I^{\lambda} f(t)$", fontsize=12)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(fontsize=11)

    fig.savefig(savefile, dpi=300)

    # 表示
    plt.show()
