# relaxation modulus of fractional Maxwell model

import numpy as np
from math import *
# from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_relaxMod(E, nu, kmax, T, t):
    relaxMod = 0
    for k in range(kmax):
        rMod = E*(-(t/T)**nu)**k / gamma(nu*k + 1)
        relaxMod += rMod

    return relaxMod

def relaxMod(E, nu, kman, T, tim):
    relaxMod = []
    for t in tim:
        r = calc_relaxMod(E, nu, kmax, T, t)
        relaxMod.append(r)
    return relaxMod

def calc_relaxFunc(nu, kman, T, t):
    relaxFunc = 0
    for k in range(kmax):
        rFunc = (-(t/T)**nu)**k / gamma(nu*k + 1)
        relaxFunc += rFunc
    relaxFunc = 1 - relaxFunc

    return relaxFunc

def relaxFunc(nu, kman, T, tim):
    relaxFunc = []
    for t in tim:
        r = calc_relaxFunc(nu, kman, T, t)
        relaxFunc.append(r)
    return relaxFunc

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
    try:
        nu = float(input('Enter fractional power (0 < nu < 1) (default = 0.5): '))
    except ValueError:
        nu = 0.5
    kappa = modulus / insMod
    relaxTime = kappa**(1/nu)*viscosity/modulus
    return insMod, modulus, relaxTime, nu

def timeAxis(relaxTime):
    log_relaxT = np.log10(relaxTime)
    linearTime = np.linspace(0, relaxTime*3, 31)
    scaledLinearTime = [t/relaxTime for t in linearTime]
    logTime = np.logspace(log_relaxT-1, log_relaxT+2, 31)
    scaledLogTime = [np.log10(t/relaxTime) for t in logTime]
    timeAxis = [linearTime, scaledLinearTime, logTime, scaledLogTime]
    return timeAxis

if __name__=='__main__':
    # calculating relaxation modulus and relaxation funcion
    insMod, modulus, relaxTime, nu = reqParams()
    param_text = '(Ei = {0:.1f} MPa, E = {1:.1f} MPa, tau = {2:.1f} ms, nu = {3:.2f})'.format(insMod/10**6, modulus/10**6, relaxTime*10**3, nu)
    timeAxis = timeAxis(relaxTime)
    try:
        select = int(input('Selection (relaxation modulus: 0, relaxation function: 1): '))
    except ValueError:
        select = 0
    
    kmax = 100

    if select == 0:
        x1 = timeAxis[0]
        x1_scaled = timeAxis[1]
        x1_label = 't/tau'
        y1 = relaxMod(insMod, nu, kmax, relaxTime, x1)
        y1_label = 'E(t) /Pa'
        label1 = 'Relaxation modulus (linear)'
        x2 = timeAxis[2]
        x2_scaled = timeAxis[3]
        x2_label = 'log[t/tau]'
        relaxMod = relaxMod(insMod, nu, kmax, relaxTime, x2)
        y2 = [np.log10(r) for r in relaxMod]
        y2_label = 'log[E(t) /Pa]'
        label2 = 'Relaxation modulus (log)'
        legend_loc='upper right'
        savefile = './png/fractional_Maxwell_relax_modulus.png'

    if select == 1:
        x1 = timeAxis[0]
        x1_scaled = timeAxis[1]
        x1_label = 't/tau'
        y1 = relaxFunc(nu, kmax, relaxTime, x1)
        y1_label = 'phi(t)'
        label1 = 'Relaxation function (linear)'
        x2 = timeAxis[2]
        x2_scaled = timeAxis[3]
        x2_label = 'log[t/tau]'
        relaxFunc = relaxFunc(nu, kmax, relaxTime, x2)
        y2 = [np.log10(r) for r in relaxFunc]
        y2_label = 'log[phi(t)]'
        label2 = 'Relaxation function (log)'
        legend_loc='upper left'
        savefile = './png/fractional_Maxwell_relax_func.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,10), tight_layout=True)
    ax1 = fig.add_subplot(211)
    ax1.set_title('fractional Maxwell model '+param_text)
    ax1.set_xlabel(x1_label)
    ax1.set_ylabel(y1_label)
#    ax1.set_xlim(-0.1,3.1)
#    ax1.set_ylim(-0.05, np.max(y1)*1.1)
    ax1.grid()
    ax1.scatter(x1_scaled, y1, c='r', label=label1)
    ax1.legend(loc=legend_loc)
    ax1.set_axisbelow(True)
    ax2 = fig.add_subplot(212)
    ax2.set_xlabel(x2_label)
    ax2.set_ylabel(y2_label)
#    ax2.set_xlim(-1.1,2.1)
#    ax2.set_ylim(np.min(y2)-1, np.max(y2)+1)
    ax2.grid()
    ax2.scatter(x2_scaled, y2, c='b', label=label2)
    ax2.legend(loc=legend_loc)
    ax2.set_axisbelow(True)
    fig.savefig(savefile)

    plt.show()