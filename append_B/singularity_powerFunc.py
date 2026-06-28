import numpy as np
import matplotlib.pyplot as plt

def calc_powerFunc(l, u):
    y = u**(l-1) * np.heaviside(u, 0.5)
    return y

if __name__=='__main__':

    l = np.linspace(-1, 2, 7)
    u = np.linspace(10**(-5), 2, 400)

    # グラフの描画
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title(r'Singularity of Power Function')
    ax.grid()
    ax.set_xlabel(r'$u$')
    ax.set_ylabel(r'$u^{{\lambda-1}}*h_s(u)$')
    for i in range(len(l)):
        y = calc_powerFunc(l[i], u)
        ax.plot(u, y, label=r'$\lambda$ = {0:.3f}'.format(l[i]))
    ax.legend(loc='upper right')
    ax.set_xlim(-0.5, 2.0)
    ax.set_ylim(0.0, 5.0)


# 表示
fig.savefig('./png/singularity_powerFunc.png', dpi=300)

plt.show()
