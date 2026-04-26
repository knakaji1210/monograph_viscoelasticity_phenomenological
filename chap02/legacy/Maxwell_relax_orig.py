# relaxation modulus of Maxwell model

import numpy as np
# from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_relaxMod(E, T, t):
    relaxMod = E*np.exp(-t/T)
    return relaxMod

def relaxMod(E, T, tim):
    relaxMod = []
    for t in tim:
        r = calc_relaxMod(E, T, t)
        relaxMod.append(r)
    return relaxMod

def calc_relaxFunc(T, t):
    relaxFunc = 1 - np.exp(-t/T)
    return relaxFunc

def relaxFunc(T, tim):
    relaxFunc = []
    for t in tim:
        r = calc_relaxFunc(T, t)
        relaxFunc.append(r)
    return relaxFunc

def reqParams():
    try:
        insMod = float(input('Enter instantaneous modulus value (MPa) (default = 1 MPa): '))*10**6
    except ValueError:
        insMod = 10**6
    try:
        viscosity = float(input('Enter viscosity value (kPa s) (default = 100 kPa s): '))*10**3
    except ValueError:
        viscosity = 10**5
    relaxTime = viscosity/insMod
    return insMod, relaxTime

def timeAxis(relaxTime):
    log_relaxT = np.log10(relaxTime)
    linearTime = np.linspace(0, relaxTime*3, 31)
    scaledLinearTime = [t/relaxTime for t in linearTime]
    logTime = np.logspace(log_relaxT-1, log_relaxT+1, 31)
    scaledLogTime = [np.log10(t/relaxTime) for t in logTime]
    timeAxis = [linearTime, scaledLinearTime, logTime, scaledLogTime]
    return timeAxis

if __name__=='__main__':
    # calculating relaxation modulus and relaxation funcion
    insMod, relaxTime = reqParams()
    param_text = '(Ei = {0:.1f} MPa, tau = {1:.1f} ms)'.format(insMod/10**6, relaxTime*10**3)
    timeAxis = timeAxis(relaxTime)
    try:
        select = int(input('Selection (relaxation modulus: 0, relaxation function: 1): '))
    except ValueError:
        select = 0
    
    if select == 0:
        x1 = timeAxis[0]
        x1_scaled = timeAxis[1]
        x1_label = 't/tau'
        y1 = relaxMod(insMod, relaxTime, x1)
        y1_label = 'E(t) /Pa'
        label1 = 'Relaxation modulus (linear)'
        x2 = timeAxis[2]
        x2_scaled = timeAxis[3]
        x2_label = 'log[t/tau]'
        relaxMod = relaxMod(insMod, relaxTime, x2)
        y2 = [np.log10(r) for r in relaxMod]
        y2_label = 'log[E(t) /Pa]'
        label2 = 'Relaxation modulus (log)'
        legend_loc='upper right'
        savefile = './png/Maxwell_relax_modulus.png'

    if select == 1:
        x1 = timeAxis[0]
        x1_scaled = timeAxis[1]
        x1_label = 't/tau'
        y1 = relaxFunc(relaxTime, x1)
        y1_label = 'phi(t)'
        label1 = 'Relaxation function (linear)'
        x2 = timeAxis[2]
        x2_scaled = timeAxis[3]
        x2_label = 'log[t/tau]'
        relaxFunc = relaxFunc(relaxTime, x2)
        y2 = [np.log10(r) for r in relaxFunc]
        y2_label = 'log[phi(t)]'
        label2 = 'Relaxation function (log)'
        legend_loc='upper left'
        savefile = './png/Maxwell_relax_func.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,10), tight_layout=True)
    ax1 = fig.add_subplot(211)
    ax1.set_title('Maxwell model '+param_text)
    ax1.set_xlabel(x1_label)
    ax1.set_ylabel(y1_label)
    ax1.set_xlim(-0.1,3.1)
    ax1.set_ylim(-np.max(y1)*0.05, np.max(y1)*1.1)
    ax1.grid()
    ax1.scatter(x1_scaled, y1, c='r', label=label1)
    ax1.legend(loc=legend_loc)
    ax1.set_axisbelow(True)
    ax2 = fig.add_subplot(212)
    ax2.set_xlabel(x2_label)
    ax2.set_ylabel(y2_label)
    ax2.set_xlim(-1.1,1.1)
    ax2.set_ylim(np.min(y2)-1, np.max(y2)+1)
    ax2.grid()
    ax2.scatter(x2_scaled, y2, c='b', label=label2)
    ax2.legend(loc=legend_loc)
    ax2.set_axisbelow(True)
    fig.savefig(savefile)

    plt.show()