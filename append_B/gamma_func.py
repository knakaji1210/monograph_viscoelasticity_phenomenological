import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp

# xの値を設定 (0より大きい範囲)
x = np.linspace(-5.0, 5.0, 400)
y = sp.gamma(x)
xd = np.array([0, 1, 2, 3, 4, 5])
yd = sp.factorial(xd, exact=True)

# グラフの描画
fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111)
ax.scatter(xd+1, yd, label=r'$(n-1)!$', color='red', zorder=1)
ax.plot(x, y, label=r'$\Gamma(x)$', color='blue', linewidth=2, zorder=0)


# グラフの装飾
ax.set_title(r'Gamma Function, $\Gamma(x)$', fontsize=14)
ax.set_xlabel('x', fontsize=12)
ax.set_ylabel(r'$\Gamma(x)$', fontsize=12)
ax.grid(True, linestyle='--', alpha=0.7)
ax.legend(loc='lower right')
ax.set_xlim(-5.0, 5.0)
ax.set_ylim(-10, 20)
ax.hlines(0, -5.0, 5.0, ls=':', color='black')
ax.vlines(0, -10.0, 20.0, ls=':', color='black')

# 表示
plt.savefig('./png/gamma_func.png', dpi=300)

plt.show()
