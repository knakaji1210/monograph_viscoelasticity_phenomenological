# Geminiとやりとりして作ったバージョンだけど結局あまりうまくいかなかった。
# 記録のために取っておく。

# relaxation modulus of fractional Maxwell model

import numpy as np
from math import *
import matplotlib.pyplot as plt
from differintP.functions import MittagLeffler
# 上のライブラリだとnuの値が小さくなるとカーブが途中で切れる
# そこでそれを使わず、MittagLeffler関数を自作する（Geminiのヘルプ付き）
# from scipy.integrate import quad
# from scipy.special import gamma

def reqParams():
    try:
        insMod = float(input('Enter instantaneous modulus value (MPa) (default = 1 MPa): '))*10**6
    except ValueError:
        insMod = 10**6
    try:
        modulus = float(input('Enter modulus value of spring-pot (MPa) (default = 0.1 MPa): '))*10**6
    except ValueError:
        modulus = 10**5
    try:
        viscosity = float(input('Enter viscosity value of spring-pot (kPa s) (default = 100 kPa s): '))*10**3
    except ValueError:
        viscosity = 10**5
    return insMod, modulus, viscosity

def timeAxes(tau):
    log_tau = np.log10(tau)
    linearTime = np.linspace(0, tau, 500)
    scaledLinearTime = linearTime/tau
    logTime = np.logspace(log_tau-1, log_tau+2, 51)
    scaledLogTime = np.log10(logTime/tau)
    timeAxes = [linearTime, scaledLinearTime, logTime, scaledLogTime]
    return timeAxes

def custom_mittag_leffler(alpha, beta, x):
    """
    数値積分と級数展開のハイブリッド実装。
    x=0付近の積分発散バグを、安定なテイラー展開で100%完全にカバーします。
    """
    # 1. 0に近い領域（-0.1 < x <= 0）は、絶対に積分を通さず級数展開で処理
    # この狭い範囲であれば、項数20で暴走（カーブの途切れ）は100%起きず、正確に1.0から減衰します。
    if abs(x) < 0.1:
        val = 0.0
        for k in range(20): 
            val += (x**k) / gamma(alpha * k + beta)
        return val
    
    # 2. x <= -0.1 の長時間領域は、分母が安定するため、滑らかな数値積分で計算
    def integrand(t):
        numerator = np.exp(-t) * np.sin(alpha * np.pi) * (t**(alpha - beta))
        denominator = (t**(2 * alpha)) + 2 * (t**alpha) * x * np.cos(alpha * np.pi) + (x**2)
        if denominator == 0:
            return 0
        return numerator / denominator
    
    # 積分範囲は元の安全な形に戻します（手前のifで0付近を弾いているため、1e-10にする必要もありません）
    integral_part, _ = quad(integrand, 0, np.inf, limit=150)
    return integral_part / np.pi

def calc_relaxMod(E, nu, x_array):
    relaxMod = E * MittagLeffler(nu, 1, x_array, num_terms=100)
    return relaxMod

def calc_relaxMod2(E, nu, x_array):
    """
    自作のミッタク＝レフラー関数を要素ごとに適用して緩和弾性率を計算
    """
    relaxMod_list = []
    for x in x_array:
        # 自作の関数を呼び出し（beta = 1 固定）
        val = E * custom_mittag_leffler(nu, 1, x)
        relaxMod_list.append(val)
        relaxMod = np.array(relaxMod_list)
    return relaxMod

cmap = plt.get_cmap('coolwarm')

if __name__=='__main__':
    # calculating relaxation modulus
    insMod, modulus, viscosity = reqParams()
    kappa = modulus / insMod
    tau = viscosity/modulus
    param_text = r'($E_i$ = {0:.1f} MPa, $E$ = {1:.1f} MPa, $\tau$ = {2:.1f} ms)'.format(insMod/10**6, modulus/10**6, tau*10**3)
    timeAxes = timeAxes(tau)
    try:
        select = int(input('Selection (relaxation modulus (linear): 0, relaxation modulus (log): 1 (default = 0): '))
    except ValueError:
        select = 0

    nu_arr = np.array([0.2, 0.4, 0.6, 0.8, 0.99])

    y_array = np.zeros((len(nu_arr), len(timeAxes[0])))  # 階数ごとの緩和スペクトルを格納する配列
    tau_prime_arr = np.zeros(len(nu_arr))  # 階数ごとの緩和時間を格納する配列

    if select == 0:
        tau_primes = [kappa**(1/nu) * tau for nu in nu_arr]
        max_tau_prime = max(tau_primes)
        tim = np.linspace(0, max_tau_prime * 5, 500)
        scaled_tim = tim / tau
        x_label = r'$t/\tau$'
        y_label = r'$E$($t$) /MPa'
        legend_loc='upper right'
        for i in range(len(nu_arr)):
            nu = nu_arr[i]
            tau_prime = tau_primes[i]
            tau_prime_arr[i] = tau_prime
            x = -(tim/tau_prime)**nu
            x[0] = 0.0
            y_array[i] = calc_relaxMod(insMod, nu, x)/10**6  # rescale to MPa
        savefile = './png/fractional_Maxwell_relaxMod_linear.png'


#        x2 = timeAxes[2]
#        x2_scaled = timeAxes[3]
#        x2_label = 'log[t/tau]'
#        relaxMod = relaxMod(insMod, nu, kmax, relaxTime, x2)
#        y2 = [np.log10(r) for r in relaxMod]
#        y2_label = 'log[E(t) /Pa]'
#        label2 = 'Relaxation modulus (log)'
#        legend_loc='upper right'
#        savefile = './png/fractional_Maxwell_relax_modulus.png'



    # drawing graphs
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title('fractional Maxwell model '+param_text)
    for i in range(len(nu_arr)):
        ax.plot(scaled_tim, y_array[i], color=cmap(i/len(nu_arr)), label='nu = {0:.2f}'.format(nu_arr[i]))
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xlim(-0.1,1.1)
    ax.set_ylim(-0.05, 1.1)
    ax.grid()
    ax.legend(loc=legend_loc, fontsize=11)
    ax.set_axisbelow(True)

#    fig.savefig(savefile, dpi=300)

    plt.show()