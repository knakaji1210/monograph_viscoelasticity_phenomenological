# マクローリン展開の計算

import math
import numpy as np
import matplotlib.pyplot as plt

def maclaurin_exp(x, n_terms):
    approx_value = 0
    for n in range(n_terms):
        approx_value += (x**n) / math.factorial(n)
    return approx_value

# xの範囲を設定
x = np.linspace(-3, 3, 200)

# 近似する次数のプロット
if __name__=='__main__':
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.grid()
    ax.plot(x, np.exp(x), label=r'Exact: $e^x$', color='black', linewidth=2)
    ax.plot(x, [maclaurin_exp(val, 4) for val in x], label=r'Maclaurin ($n$ = 3)', linestyle='--')
    ax.plot(x, [maclaurin_exp(val, 6) for val in x], label=r'Maclaurin ($n$ = 5)', linestyle='--')
    ax.plot(x, [maclaurin_exp(val, 11) for val in x], label=r'Maclaurin ($n$ = 10)', linestyle='--')

    ax.set_ylim(-2, 10)
    ax.set_xlim(-3, 3)
    ax.legend()
    ax.set_title(r'Maclaurin Expansion of $e^x$')
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.grid(True)

    fig.savefig('./png/MaclaurinExpansion.png', dpi=300)

    plt.show()
