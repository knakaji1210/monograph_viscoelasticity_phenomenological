# relaxation modulus of fractional Maxwell model

import numpy as np
from math import *
import matplotlib.pyplot as plt
from differintP.functions import MittagLeffler

def reqParams():
    try:
        insMod = float(input('Enter instantaneous modulus value (MPa) (default = 1 MPa): '))*10**6
    except ValueError:
        insMod = 10**6
    try:
        modulus = float(input('Enter modulus value of spring-pot (MPa) (default = 0.8 MPa): '))*10**6
    except ValueError:
        modulus = 8*10**5
    try:
        viscosity = float(input('Enter viscosity value of spring-pot (kPa s) (default = 100 kPa s): '))*10**3
    except ValueError:
        viscosity = 10**5
    return insMod, modulus, viscosity

def timeAxes(tau):
    log_tau = np.log10(tau)
    linearTime = np.linspace(10**(-5), tau*3, 500)
    scaledLinearTime = linearTime/tau
    logTime = np.logspace(log_tau-2, log_tau+0.5, 500)
    scaledLogTime = np.log10(logTime/tau)
    timeAxes = [linearTime, scaledLinearTime, logTime, scaledLogTime]
    return timeAxes

def calc_relaxMod(E, nu, x):
    # E = Ei
    relaxMod = E*MittagLeffler(nu, 1, x, num_terms=300)
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

    nu_arr = np.array([0.2, 0.4, 0.6, 0.8, 1.0])

    y_array = np.zeros((len(nu_arr), len(timeAxes[0])))  # 階数ごとの緩和スペクトルを格納する配列
    tau_prime_arr = np.zeros(len(nu_arr))  # 階数ごとの緩和時間を格納する配列

    if select == 0:
        tim = timeAxes[0]
        scaled_tim = timeAxes[1]
        x_label = r'$t/\tau$'
        y_label = r'$E$($t$) /MPa'
        for i in range(len(nu_arr)):
            nu = nu_arr[i]
            tau_prime = kappa**(1/nu)*tau
            tau_prime_arr[i] = tau_prime
            x = -(tim/tau_prime)**nu
            y_array[i] = calc_relaxMod(insMod, nu, x)/10**6  # rescale to MPa
        x_lim = [-0.1, scaled_tim[-1]*1.05]
        y_lim = [-0.05*np.max(y_array), 1.1*np.max(y_array)]
        legend_loc='upper right'
        savefile = './png/fractional_Maxwell_relax_modulus_linear.png'

    if select == 1:
        tim = timeAxes[2]
        scaled_tim = timeAxes[3]
        x_label = r'log($t/\tau$)'
        y_label = r'log($E(t)$ /Pa)'
        for i in range(len(nu_arr)):
            nu = nu_arr[i]
            tau_prime = kappa**(1/nu)*tau
            tau_prime_arr[i] = tau_prime
            x = -(tim/tau_prime)**nu
            y_array[i] = np.log10(calc_relaxMod(insMod, nu, x))
        x_lim = [scaled_tim[0], scaled_tim[-1]]
        y_lim = [np.min(y_array), 1.05*np.max(y_array)]
        legend_loc='upper right'
        savefile = './png/fractional_Maxwell_relax_modulus_log.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title('fractional Maxwell model '+param_text)
    for i in range(len(nu_arr)):
        ax.plot(scaled_tim, y_array[i], color=cmap(i/len(nu_arr)), label=r'$\nu$ = {0:.1f}'.format(nu_arr[i]))
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xlim(x_lim[0], x_lim[1])
    ax.set_ylim(y_lim[0], y_lim[1])
    ax.grid()
    ax.legend(loc='upper right', fontsize=11)
    ax.set_axisbelow(True)

    fig.savefig(savefile, dpi=300)

    plt.show()