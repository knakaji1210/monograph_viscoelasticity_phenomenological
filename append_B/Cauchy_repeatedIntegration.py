# 「反復積分に関するコーシーの公式」を物理的にイメージするための試み

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp

cmap = plt.get_cmap('coolwarm')

def calc_originalFunc(u):
    y = -(u-1)**2 + 2*(u-1)
    return y

def calc_filterFunc(u, t):
    y = t-u
    return y

if __name__=='__main__':

    t = np.array([1.2, 2.0, 2.8])
    u = np.linspace(0, 3.0, 400)

    # グラフの描画
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title('Cauchy Repeated Integration', fontsize=14)
    ax.grid()
    ax.set_xlabel(r'$u$')
    ax.set_ylabel(r'$f(u)$, $(t-u)$, $(t-u)f(u)$')
    ax.plot(u, calc_originalFunc(u), label=r'$f(u)$', color='blue', linewidth=2)
    for i in range(len(t)):
        y1 = calc_filterFunc(u, t[i])
        y2 = calc_originalFunc(u)*calc_filterFunc(u, t[i])
        ax.plot(u, y1, label=r'$(t-u)$, $t$ = {0:.1f}'.format(t[i]), color=cmap((i+1)*0.2), linewidth=1, ls='--')
        ax.plot(u, y2, label=r'$(t-u)f(u)$, $t$ = {0:.1f}'.format(t[i]), color=cmap((i+1)*0.2), linewidth=1, ls='-')
    ax.legend(loc='upper right')
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 1.25)


# 表示
fig.savefig('./png/Cauchy_repeatedIntegration.png', dpi=300)

plt.show()
