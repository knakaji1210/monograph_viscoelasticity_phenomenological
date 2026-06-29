import numpy as np
import matplotlib.pyplot as plt
from differintP import CaputoL1point  # differintの場合も同様

# 1. 時間領域と関数の定義
t_start = 0.0
t_end = 5.0
N = 100
t_values = np.linspace(t_start, t_end, N)

# 対象とする関数（例: f(t) = t^2）
def my_func(t):
    return np.sin(t)

# 2. Caputo微分の計算
alpha = 0.5
caputo_derivs = []

for t in t_values:
    if t == t_start:
        val = 0.0
    else:
        # 引数の順番を修正: (alpha, 関数, 開始点, 終了点, 分割数)
        val = CaputoL1point(alpha, my_func, t_start, t, N)
    caputo_derivs.append(val)

# 3. グラフのプロット
plt.figure(figsize=(8, 5))
plt.plot(t_values, my_func(t_values), label='Original: f(t) = t^2', color='blue')
plt.plot(t_values, caputo_derivs, label=f'Caputo Derivative (alpha={alpha})', color='red', linestyle='--')
plt.title('Caputo Fractional Derivative using differintP')
plt.xlabel('t')
plt.ylabel('Value')
plt.legend()
plt.grid(True)
plt.show()
