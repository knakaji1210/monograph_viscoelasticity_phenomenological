# Levy Stable Distributionの確率密度関数(PDF)の描画

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import levy_stable

# 描画するx軸の範囲
x = np.linspace(-5, 5, 1000)

# パラメータパターンをリストで定義 (alpha, beta)
params = [
    (2.0, 0.0),  # 正規分布
    (1.0, 0.0),  # コーシー分布
    (1.5, 0.0)   # 対称レヴィ分布
]

func_list = [r'Gaussian ($\alpha=2.0, \beta=0.0$)', r'Cauchy ($\alpha=1.0, \beta=0.0$)', r'Levy ($\alpha=1.5, \beta=0.0$)']  # 関数のリスト

# グラフの設定
fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111)

# 各パラメータで確率密度関数(PDF)を計算してプロット
for alpha, beta in params:
    # levy_stable.pdf(x, alpha, beta, loc=位置母数, scale=尺度母数)
    pdf = levy_stable.pdf(x, alpha, beta)
    ax.plot(x, pdf, label=func_list[params.index((alpha, beta))])

# グラフの装飾
ax.set_title('Probability Density Function (PDF)', fontsize=14)
ax.set_xlabel('$x$', fontsize=12)
ax.set_ylabel(f'$\\rho(x)$', fontsize=12)
ax.set_xlim(-5, 5)
ax.set_ylim(0, 0.4)
ax.legend(fontsize=12)
ax.grid(True, linestyle='--', alpha=0.7)

savefile = './png/Levy_stable_distribution2.png'
# グラフの保存
fig.savefig(savefile, dpi=300)

# グラフの表示
plt.show()
