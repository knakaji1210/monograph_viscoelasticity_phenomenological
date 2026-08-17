import numpy as np
import matplotlib.pyplot as plt

def reqParams():
    try:
        insMod = float(input('Enter instantaneous modulus value (GPa) (default = 1 GPa): '))*10**9
    except ValueError:
        insMod = 10**9
    try:
        infMod = float(input('Enter equilibrium modulus value (MPa) (default = 1 MPa): '))*10**6
    except ValueError:
        infMod = 10**6
    try:
        modulus = float(input('Enter modulus value of spring-pot (GPa) (default = 0.1 GPa): '))*10**9
    except ValueError:
        modulus = 10**8
    try:
        viscosity = float(input('Enter viscosity value of spring-pot (kPa s) (default = 1000 kPa s): '))*10**3
    except ValueError:
        viscosity = 10**6
    return insMod, infMod, modulus, viscosity

def relaxFunc(x, nu, k):
    kai = np.cos(np.pi*nu/2)
    return nu*(1-1/k)*x*(1+2*x+kai*x**2)/(1+2*kai*x+x**2)**2

cmap = plt.get_cmap('coolwarm')

if __name__=='__main__':
    # calcul1ating relaxation spectrum
    insMod, infMod, modulus, viscosity = reqParams()
    k = insMod/infMod
    kappa = modulus / (insMod - infMod)
    tau = viscosity/modulus
    param_text = """
    ($E_i$ = {0:.1f} MPa, $E_{{\infty}}$ = {1:.1f} MPa, $E$ = {2:.1f} MPa, $\\tau$ = {3:.1f} ms)""".format(insMod/10**6, infMod/10**6, modulus/10**6, tau*10**3)
    tim = np.logspace(int(np.log10(tau))-10, int(np.log10(tau))+4, 500)
    tim_label = r'log($\tau$ / s)'    # 横軸をtimではなくtauと書きたいからこう書いている
    log_Tim = np.log10(tim)
    ylabel = r'$H(\tau) / E_i$'

    nu_arr = np.array([0.2, 0.4, 0.6, 0.8, 1.0])

    y_array = np.zeros((len(nu_arr), len(tim)))  # 階数ごとの緩和スペクトルを格納する配列
    tau_prime_arr = np.zeros(len(nu_arr))  # 階数ごとの緩和時間を格納する配列

    for i in range(len(nu_arr)):
        nu = nu_arr[i]
        tau_prime = kappa**(1/nu)*tau
        tau_prime_arr[i] = tau_prime
        x = kappa*(tau/tim)**nu
        y_array[i] = relaxFunc(x, nu, k)

    fig = plt.figure(tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title('Relaxation spectrum of fractional Zener model '+param_text)
    ax.set_ylim(-0.1,1.1)
    ax.set_xlabel(tim_label)
    ax.set_ylabel(ylabel)
    for i in range(len(nu_arr)):
        ax.plot(log_Tim, y_array[i], label=f'$\\nu$ = {nu_arr[i]:.1f}', color=cmap(i/len(nu_arr)), linewidth=2, zorder=2)
        ax.vlines([np.log10(tau_prime_arr[i])], 0, 1.1, label=f'$\\tau^\\prime$ = {tau_prime_arr[i]*10**6:.1f} $\\mu$s', color=cmap(i/len(nu_arr)), ls='--', linewidth=1, zorder=1)
    ax.vlines([np.log10(tau)], 0, 1.1, label=f'$\\tau$ = {tau*10**3:.1f} ms', color='red', ls='--', linewidth=1, zorder=1)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(loc='upper left', fontsize=11)

    savefile = './png/Alfrey_approx_fracZener_relaxSpectrum_strMod.png'
    fig.savefig(savefile, dpi=300)

    plt.show()