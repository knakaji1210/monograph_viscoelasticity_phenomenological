import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp

# xの値を設定 (0より大きい範囲)
x = np.linspace(-5.0, 5.0, 400)
y = sp.gamma(x)
xd = np.array([0, 1, 2, 3, 4, 5])
yd = sp.factorial(xd, exact=True)

# グラフの描画
plt.figure(figsize=(8, 6))
plt.scatter(xd+1, yd, label=r'factrial($n$)', color='red', zorder=1)
plt.plot(x, y, label=r'$\Gamma(x)$', color='blue', linewidth=2, zorder=0)


# グラフの装飾
plt.title('Gamma Function', fontsize=14)
plt.xlabel('x', fontsize=12)
plt.ylabel(r'$\Gamma(x)$', fontsize=12)
plt.axhline(0, color='black', linewidth=1)
plt.axvline(0, color='black', linewidth=1)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=12)
plt.ylim(-10, 20)

# 表示
plt.savefig('./png/gamma_func.png', dpi=300)

plt.show()
