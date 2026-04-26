# creep compliance of Kelvin-Voigt model

import numpy as np
# from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_creepComp(E, T, t):
    creepComp = (1 - np.exp(-t/T))/E
    return creepComp

def creepComp(E, T, tim):
    creepComp = []
    for t in tim:
        r = calc_creepComp(E, T, t)
        creepComp.append(r)
    creepComp = np.array(creepComp)
    return creepComp

def calc_creepFunc(T, t):
    creepFunc = 1 - np.exp(-t/T)
    return creepFunc

def creepFunc(T, tim):
    creepFunc = []
    for t in tim:
        r = calc_creepFunc(T, t)
        creepFunc.append(r)
    creepFunc = np.array(creepFunc)
    return creepFunc

def reqParams():
    try:
        infMod = float(input('Enter equilibrium modulus value (MPa) (default = 1 MPa): '))*10**6
    except ValueError:
        infMod = 10**6
    try:
        viscosity = float(input('Enter viscosity value (kPa s) (default = 100 kPa s): '))*10**5
    except ValueError:
        viscosity = 10**5
    retardTime = viscosity/infMod
    return infMod, retardTime

def timeAxis(retardTime):
    log_relaxT = np.log10(retardTime)
    linearTime = np.linspace(0, retardTime*3, 31)
    scaledLinearTime = [t/retardTime for t in linearTime]
    logTime = np.logspace(log_relaxT-1, log_relaxT+1, 31)
    scaledLogTime = [np.log10(t/retardTime) for t in logTime]
    timeAxis = [linearTime, scaledLinearTime, logTime, scaledLogTime]
    return timeAxis

if __name__=='__main__':
    # calculating creep compliance and creep funcion
    infMod, retardTime = reqParams()
    param_text = r'($E_{{inf}}$ = {0:.1f} MPa, $\tau$ = {1:.1f} ms)'.format(infMod/10**6, retardTime*10**3)
    timeAxis = timeAxis(retardTime)
    try:
        select = int(input('Selection (creep compliance: 0, creep function: 1): '))
    except ValueError:
        select = 0
    
    if select == 0:
        x1 = timeAxis[0]
        x1_scaled = timeAxis[1]
        x1_label = r'$t$/$\tau$'
        y1 = creepComp(infMod, retardTime, x1)
        y1 *= 10**6             # rescale to MPa
        y1_label = r'$J$($t$) /MPa$^{{-1}}$'
        label1 = 'Creep compliance (linear)'
        x2 = timeAxis[2]
        x2_scaled = timeAxis[3]
        x2_label = r'log[$t$/$\tau$]'
        creepComp = creepComp(infMod, retardTime, x2)
#       y2 = [np.log10(r) for r in creepComp]   # 初期のバージョン
        y2 = np.log10(creepComp)                 # creepCompの出力をnp.ndarrayに変更したのでこのように簡単に書ける
        y2_label = r'log[$J$($t$) /Pa$^{{-1}}$]'
        label2 = 'Creep compliance (log)'
        legend_loc='upper left'
        savefile = './png/Voigt_creep_compliance.png'

    if select == 1:
        x1 = timeAxis[0]
        x1_scaled = timeAxis[1]
        x1_label = r'$t$/$\tau$'
        y1 = creepFunc(retardTime, x1)
        y1_label = r'$\psi$($t$)'
        label1 = 'Creep function (linear)'
        x2 = timeAxis[2]
        x2_scaled = timeAxis[3]
        x2_label = r'log[$t$/$\tau$]'
        creepFunc = creepFunc(retardTime, x2)
#       y2 = [np.log10(r) for r in creepFunc]   # 初期のバージョン
        y2 = np.log10(creepFunc)                 # creepFuncの出力をnp.ndarrayに変更したのでこのように簡単に書ける
        y2_label = r'log[$\psi$($t$)]'
        label2 = 'Creep function (log)'
        legend_loc='upper left'
        savefile = './png/Voigt_creep_func.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,10), tight_layout=True)
    ax1 = fig.add_subplot(211)
    ax1.set_title('Kelvin-Voight model '+param_text)
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
    fig.savefig(savefile, dpi=300)

    plt.show()