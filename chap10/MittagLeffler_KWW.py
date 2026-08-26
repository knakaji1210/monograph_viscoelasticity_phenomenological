# relaxation modulus of fractional Maxwell model

import numpy as np
from math import *
import matplotlib.pyplot as plt
from differintP.functions import MittagLeffler
import scipy.special as sp

def reqParams():
    try:
        insMod = float(input('Enter instantaneous modulus value (MPa) (default = 0.1 GPa): '))*10**9
    except ValueError:
        insMod = 10**8
    try:
        infMod = float(input('Enter equilibrium modulus value (MPa) (default = 1 MPa): '))*10**6
    except ValueError:
        infMod = 10**6
    try:
        modulus = float(input('Enter modulus value of spring-pot (MPa) (default = 0.1 MPa): '))*10**6
    except ValueError:
        modulus = 10**7
    try:
        viscosity = float(input('Enter viscosity value of spring-pot (kPa s) (default = 1000 kPa s): '))*10**3
    except ValueError:
        viscosity = 10**6
    try:
        nu = float(input('Enter fractional power (0 < nu < 1) (default = 0.5): '))
    except ValueError:
        nu = 0.5
    return insMod, infMod, modulus, viscosity, nu

def timeAxes(tau):
    # tau must be tau_prime
    linearTime = np.linspace(10**(-8), tau*2, 500)
    scaledLinearTime = linearTime/tau
    timeAxes = [linearTime, scaledLinearTime]
    return timeAxes

def calc_creepComp_ML(k, nu, x):
    creepComp_ML = 1 - (1 - 1/k) * MittagLeffler(nu, 1, x, num_terms=300)
    return creepComp_ML

def calc_creepComp_KWW(k, x):
    creepComp_KWW = 1 - (1 - 1/k) * np.exp(x)
    return creepComp_KWW

if __name__=='__main__':
    # calculating creep compliance of fractional Zener model
    insMod, infMod, modulus, viscosity, nu = reqParams()
    kappa = modulus / (insMod - infMod)
    k = insMod / infMod
    print('kappa = {0:.3f}, k = {1:.3f}'.format(kappa, k))
    tau = viscosity/modulus
    print(kappa**(1/nu))
    tau_prime = kappa**(1/nu)*tau
    tau1 = tau_prime * sp.gamma(1+nu)**(1/nu)
    print('tau = {0:.3f} s, tau\' = {1:.10f} s, tau1 = {2:.10f} s'.format(tau, tau_prime, tau1))
    param_text = r'($E_i$ = {0:.1f} MPa, $E_\infty$ = {1:.1f} MPa, $E$ = {2:.1f} MPa, $\tau$ = {3:.2f} s)'.format(insMod/10**6, infMod/10**6, modulus/10**6, tau)
    timeAxes = timeAxes(tau_prime)
    tim = timeAxes[0]
    scaled_tim = timeAxes[1]
    x_label = r'$t/\tau^{{\prime}}$'
    y_label = r'$\phi$($t$)'
    x_ML = -(tim/tau_prime)**nu
    y_ML = calc_creepComp_ML(k, nu, x_ML)
    x_KWW = -(tim/tau1)**nu
    y_KWW = calc_creepComp_KWW(k, x_KWW)
    x_SLS = -(tim/tau_prime)
    y_SLS = calc_creepComp_KWW(k, x_SLS)

    x_lim = [-0.1, scaled_tim[-1]*1.05]
    y_lim = [0,1]
    legend_loc='upper right'
    savefile = './png/MittagLeffler_KWW.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title('fractional Zener model '+param_text)
    ax.plot(scaled_tim, y_ML, label=r'Mittag-Leffler ($\nu$ = {0:.1f}, $\tau^{{\prime}}$ = {1:.2f} ms)'.format(nu, tau_prime*10**3))
    ax.plot(scaled_tim, y_KWW, label=r'KWW ($\nu$ = {0:.1f}, $\tau_1$ = {1:.2f} ms)'.format(nu, tau1*10**3))
    ax.plot(scaled_tim, y_SLS, label=r'SLS ($\nu$ = 1.0, $\tau^{{\prime}}$ = {0:.2f} ms)'.format(tau_prime*10**3))
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xlim(x_lim[0], x_lim[1])
    ax.set_ylim(y_lim[0], y_lim[1])
    ax.grid()
    ax.legend(loc='upper left', fontsize=11)
    ax.set_axisbelow(True)

    fig.savefig(savefile, dpi=300)

    plt.show()