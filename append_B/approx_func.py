import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp

def calc_approxFunc(l, t):
    y = t**(l-1) * np.heaviside(t, 0.5) / sp.gamma(l)
    return y

if __name__=='__main__':

    l = np.array([0.001, 0.01, 0.1])
    t = np.linspace(10**(-5), 1.0, 400)

    # グラフの描画
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.grid()
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$t^{{\lambda-1}}$/$\Gamma(\lambda)$')
    for i in range(len(l)):
        y = calc_approxFunc(l[i], t)
        ax.plot(t, y, label=r'$\lambda$ = {0:.3f}'.format(l[i]))
    ax.legend(loc='upper right')
    ax.set_ylim(0, 5)


# 表示
fig.savefig('./png/approx_func.png', dpi=300)

plt.show()
