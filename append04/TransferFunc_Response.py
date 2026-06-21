import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# パラメータの設定
tau = 0.1       # 時定数
E = 1.0         # 入力の大きさ（ここでは弾性率/MPa）
k = 10

select_text = 'Selection (Maxwell: 0, SLS2: 1): '  
try:
    select = int(input(select_text))
except ValueError:
    select = 0

# 1. 伝達関数の定義
# 分子 (num) と分母 (den) の係数を降べきの順で指定
if select == 0: # Maxwell model
    # G(s) = E*s / (s + 1/tau)
    num = [E, 0]            # E, 0 は s の係数と定数項
    den = [1, 1/tau]        # 1, 1/tau は s の係数と定数項
    model_name = 'Maxwell'
    param_text = r'($E_i$ = {0:.1f} MPa, $\tau$ = {1:.1f} ms)'.format(E, tau*10**3)
    savefile = './png/Maxwell_TransferFunc.png'
elif select == 1: # SLS2 model
    # G(s) = E*(s + 1/tau*k) / (s + 1/tau)
    num = [E, E/(tau*k)]    # E, E/(tau*k) は s の係数と定数項
    den = [1, 1/tau]        # 1, 1/tau は s の係数と定数項
    model_name = 'SLS II'
    param_text = r'($E_i$ = {0:.1f} MPa, $\tau$ = {1:.1f} ms, $k$ = {2:.1f})'.format(E, tau*10**3, k)
    savefile = './png/SLSII_TransferFunc.png'
sys = signal.TransferFunction(num, den)

# 2. 時間軸の作成 (0秒から1秒まで100分割)
t = np.linspace(0, 1, 100)

# 3. ステップ応答の計算
t, y = signal.step(sys, T=t)

# 4. グラフ描画
fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111)
ax.set_title(model_name + ' model ' + param_text)
ax.set_xlabel(r'$t$ /s')
ax.set_ylabel(r'$E(t)$ /MPa')
ax.grid(True)
ax.plot(t, y, label='Response to step input', color='blue')
ax.legend(loc='upper right')

plt.savefig(savefile, dpi=300)

plt.show()
