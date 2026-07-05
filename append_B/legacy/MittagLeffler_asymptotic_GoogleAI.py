import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma

def mittag_leffler_safe(alpha, beta, x, num_terms=100):
    """
    小さい |x| には級数展開、大きい負の x には漸近展開を切り替えて計算する
    """
    x = np.asarray(x, dtype=float)
    result = np.zeros_like(x)
    
    threshold = 5.0
    
    # 1. 小さい領域：級数展開
    idx_small = np.abs(x) <= threshold
    if np.any(idx_small):
        x_s = x[idx_small]
        res_s = np.zeros_like(x_s)
        for k in range(num_terms):
            log_term = k * np.log(np.abs(x_s) + 1e-15) - np.log(gamma(alpha * k + beta))
            term = np.exp(log_term) * np.sign(x_s)**k
            res_s += term
        result[idx_small] = res_s
        
    # 2. 大きい負の領域：漸近展開
    idx_large = (x < -threshold)
    if np.any(idx_large):
        x_l = x[idx_large]
        res_l = np.zeros_like(x_l)
        for k in range(1, 20):
            g_val = gamma(beta - alpha * k)
            if np.isinf(g_val) or np.isnan(g_val):
                continue
            term = 1.0 / (g_val * (x_l**k))
            res_l -= term
        result[idx_large] = res_l
        
    return result

if __name__ == '__main__':
    # パラメータ設定
    a = -1.0
    nu = 0.7  # 伸縮指数
    tau = 1.0 # 緩和時間

    # 時間軸の設定（10^-2 から 10^4 まで）
    log10_ts = -2
    log10_te = 4
    n = 1000
    t_vals = np.logspace(log10_ts, log10_te, n)

    # 1. Mittag-Leffler関数の計算
    x_vals = a * (t_vals / tau)**nu
    y_ml = mittag_leffler_safe(nu, 1.0, x_vals)
    
    # 2. KWW関数の計算
    y_kww = np.exp(-(t_vals / tau)**nu)
    
    # 3. 長時間側の極限（べき関数）の計算
    # t が大きい領域でのみ意味を持つため、プロット時に形状が分かりやすいよう全域で計算
    y_asymp = 1.0 / (gamma(1.0 - nu) * (t_vals / tau)**nu)

    # グラフの描画
    fig, ax = plt.subplots(figsize=(8, 5), tight_layout=True)
    
    # 各関数のプロット
    ax.plot(t_vals, y_ml, label=f'Mittag-Leffler ($E_{{{nu}}}(-t^{{{nu}}})$)', color='blue', linewidth=2.5)
    ax.plot(t_vals, y_kww, label=f'KWW ($\exp(-t^{{{nu}}})$)', color='orange', linestyle='--')
    ax.plot(t_vals, y_asymp, label=f'Long-time Limit ($\\propto t^{{-{nu}}}$)', color='crimson', linestyle=':', linewidth=2)
    
    # 軸の設定
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('t', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    
    # グラフの見栄え調整（y軸の下限を少し下げてべき関数の直線性を見やすくします）
    ax.set_ylim(1e-5, 2) 
    
    ax.grid(True, which="both", linestyle=":", alpha=0.6)
    ax.legend(fontsize=11)

    plt.show()
