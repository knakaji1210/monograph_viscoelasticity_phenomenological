import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma

def caputo_l1_derivative(t, f_values, alpha):
    """Caputo分数階微分 (L1算法)"""
    n = len(t)
    dt = t[1] - t[0] # 刻み幅
    caputo_df = np.zeros(n)
    caputo_df[0] = 0.0 # 初期値は0
    
    for i in range(1, n):
        s = 0.0
        for j in range(i):
            df = f_values[j + 1] - f_values[j]
            weight = (i - j) ** (1 - alpha) - (i - j - 1) ** (1 - alpha)
            s += weight * df
        caputo_df[i] = s / (gamma(2 - alpha) * (dt ** alpha))
    return caputo_df

def rl_gl_derivative(t, f_values, alpha):
    """リーマン・リウヴィル分数階微分 (GL算法)"""
    n = len(t)
    dt = t[1] - t[0] # 刻み幅
    rl_df = np.zeros(n)
    
    # 重み係数の計算
    w = np.zeros(n)
    w[0] = 1.0
    for j in range(1, n):
        w[j] = w[j-1] * (1 - (alpha + 1) / j)
        
    for i in range(n):
        s = 0.0
        for j in range(i + 1):
            s += w[j] * f_values[i - j]
        rl_df[i] = s / (dt ** alpha)
    return rl_df

# 1. データの定義
t_start = 0
t_end = 2
num_points = 200 # 始点付近を滑らかにするため少し細かくします
t_vals = np.linspace(t_start, t_end, num_points)

# 関数の定義：違いが分かりやすいように初期値を0.5にずらします
# f(t) = t^2 + 0.5
f_values = np.sin(t_vals*4) + 0.5

# 2. 0.5階微分の計算
alpha = 0.5
caputo_res = caputo_l1_derivative(t_vals, f_values, alpha)
rl_res = rl_gl_derivative(t_vals, f_values, alpha)

# 3. グラフプロット
plt.figure(figsize=(10, 6))

# 元の関数
plt.plot(t_vals, f_values, label="Original: f(t) = t^2 + 0.5", color="black", linestyle="--", alpha=0.7)

# Caputo微分
plt.plot(t_vals, caputo_res, label=f"Caputo (alpha={alpha})", color="blue", linewidth=2)

# Riemann-Liouville微分 (t=0付近の極端な発散を見やすくするためy軸を制限します)
plt.plot(t_vals, rl_res, label=f"Riemann-Liouville (alpha={alpha})", color="red", linewidth=2)

plt.xlim(-0.1, 2.1)
plt.ylim(-0.5, 5.0) # 見やすさのためにy軸の表示範囲を設定
plt.xlabel("t")
plt.ylabel("y")
plt.title("Comparison: Caputo vs Riemann-Liouville Fractional Derivative")
plt.legend()
plt.grid(True)
plt.show()
