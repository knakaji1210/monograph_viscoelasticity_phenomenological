import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# 1. 記号と関数の定義
t, s = sp.symbols('t s', real=True, positive=True)
f_t = sp.sin(t)  # ここに変換したい関数を定義

# 2. ラプラス変換の実行
# laplace_transformは (変換後の式, 収束域, 条件) を返すため、[0]で式のみを取得
F_s = sp.laplace_transform(f_t, t, s)[0]

print(f"元の関数 f(t) = {f_t}")
print(f"ラプラス変換 F(s) = {F_s}")

# 3. SymPyの数式をNumPyで計算できる関数に変換 (lambdify)
f_num = sp.lambdify(t, f_t, 'numpy')
F_num = sp.lambdify(s, F_s, 'numpy')

# 4. グラフ描画用のデータ生成
t_vals = np.linspace(0, 10, 400)
s_vals = np.linspace(0.1, 5, 400)  # s=0での発散を避けるため0.1から開始

f_y = f_num(t_vals)
F_y = F_num(s_vals)

# 5. グラフのプロット
plt.figure(figsize=(12, 5))

# 左側：元の関数 f(t)
plt.subplot(1, 2, 1)
plt.plot(t_vals, f_y, label=f'$f(t) = {sp.latex(f_t)}$', color='blue')
plt.title('Original Function $f(t)$')
plt.xlabel('t')
plt.ylabel('f(t)')
plt.grid(True)
plt.legend()

# 右側：ラプラス変換後の関数 F(s)
plt.subplot(1, 2, 2)
plt.plot(s_vals, F_y, label=f'$F(s) = {sp.latex(F_s)}$', color='red')
plt.title('Laplace Transform $F(s)$')
plt.xlabel('s')
plt.ylabel('F(s)')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()
