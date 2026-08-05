import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt

def calc_powerFunc(l, t):
    y = t**(l-1) / sp.gamma(l) * np.heaviside(t, 0.5)
    return y

if __name__=='__main__':

    l = np.linspace(-1, 2, 7)  # λの値を設定
    t = np.linspace(10**(-5), 2, 400)

    # グラフの描画
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title(r'Kernel Function of Riemann-Liouville Calculus', fontsize=14)
    ax.grid()
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$\Phi_{\lambda}(t) = t^{{\lambda-1}}/\Gamma(\lambda)*h_s(t)$')
    for i in range(len(l)):
        y = calc_powerFunc(l[i], t)
        ax.plot(t, y, label=r'$\lambda$ = {0:.3f}'.format(l[i]))
    ax.legend(loc='upper right')
    ax.set_xlim(-0.5, 2.0)
    ax.set_ylim(-5.0, 5.0)


# 表示
fig.savefig('./png/kernelFunc_RL.png', dpi=300)

plt.show()
