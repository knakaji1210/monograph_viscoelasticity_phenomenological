import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp

def calc_approxFunc(e, u):
    y = u**(e-1) * np.heaviside(u, 0.5) / sp.gamma(e)
    return y

if __name__=='__main__':

    e = np.array([0.001, 0.01, 0.1])
    u = np.linspace(10**(-5), 1.0 - 10**(-5), 400)

    # グラフの描画
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title(r'Approximated Delta Function')
    ax.grid()
    ax.set_xlabel(r'$u$')
    ax.set_ylabel(r'$u^{{\epsilon-1}}$/$\Gamma(\epsilon)$')
    for i in range(len(e)):
        y = calc_approxFunc(e[i], u)
        ax.plot(u, y, label=r'$\epsilon$ = {0:.3f}'.format(e[i]))
    ax.legend(loc='upper right')
    ax.set_xlim(-0.2, 1)
    ax.set_ylim(0, 5)


# 表示
fig.savefig('./png/approx_deltaFunc.png', dpi=300)

plt.show()
