# creep compliance of SLS I model

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def creepComp(E, k, tau, t_elapsed):
    creepComp = np.where(t_elapsed >= 0, (1 - (1 - 1/k)*np.exp(-t_elapsed/tau))/E, 0)
    # E must be infMod
    return creepComp

def creepFunc(k, tau, t_elapsed):
    creepFunc = np.where(t_elapsed >= 0, -1 + k*(1 - (1 - 1/k)*np.exp(-t_elapsed/tau)), 0)
    #creepFunc = k*(1 - (1 - 1/k)*np.exp(-t_elapsed/tau))   #Johnsonの定義
    return creepFunc

def reqParams():
    try:
        insMod = float(input('Enter instantaneous modulus value (MPa) (default = 1.0 MPa): '))*10**6
    except ValueError:
        insMod = 10**6
    try:
        infMod = float(input('Enter equilibrium modulus value (MPa) (default = 0.5 MPa): '))*10**6
    except ValueError:
        infMod = 5*10**5
    try:
        viscosity = float(input('Enter viscosity value (kPa s) (default = 500.0 kPa s): '))*10**3
    except ValueError:
        viscosity = 5*10**5
    modulus = insMod*infMod/(insMod - infMod)
    retardTime = viscosity/modulus
    k = insMod/infMod
    return insMod, infMod, k, retardTime

def timeAxis(retardTime):
    log_relaxT = np.log10(retardTime)
    linearTime = np.linspace(-retardTime*0.5, retardTime*3, 400)
    scaledLinearTime = linearTime/retardTime
    logTime = np.logspace(log_relaxT-2.0, log_relaxT+2.0, 400)
    scaledLogTime = np.log10(logTime/retardTime)
    return linearTime, scaledLinearTime, logTime, scaledLogTime

def fitTimes():
    try:
        minTime = float(input('Enter minimum time for fitting in log scale (default = -0.5): '))
    except ValueError:
        minTime = -0.5
    try:
        maxTime = float(input('Enter maximum time for fitting in log scale (default = 0.2): '))
    except ValueError:
         maxTime = 0.2
    fitTimes = [minTime, maxTime]
    return fitTimes

def getNearestIdx(list, num):
    idx = np.abs(np.asarray(list) - num).argmin()
    return idx

def fitRegion(time, minNum, maxNum):
    minFit = getNearestIdx(time, minNum)
    maxFit = getNearestIdx(time, maxNum)
    fitRegion = [minFit, maxFit]
    return fitRegion

def loglogFit(x, a, b):
    return  a*x + b

def fittedArray(x_array, param):
    fitted_array = [loglogFit(num, param[0], param[1]) for num in x_array]
    return fitted_array

def curveFit(x, y, fitTimes):
        minFit = fitRegion(x, fitTimes[0], fitTimes[1])[0]
        maxFit = fitRegion(x, fitTimes[0], fitTimes[1])[1]
        param,_ = curve_fit(loglogFit, x[minFit:maxFit], y[minFit:maxFit])
        y_fit = fittedArray(x, param)
        return y_fit, param

if __name__=='__main__':
    # calculating creep compliance and creep funcion
    insMod, infMod, k, retardTime = reqParams()
    param_text = r'($E_i$ = {0:.1f} MPa, $E_\infty$ = {1:.1f} MPa, $\tau$ = {1:.2f} s)'.format(insMod/10**6, infMod/10**6, retardTime)
    timeAxis = timeAxis(retardTime)
    try:
        select = int(input('Selection (creep compliance: 0, creep function: 1): '))
    except ValueError:
        select = 0

    try:
        fitting = int(input('Selection (curve fit: 0, no curve fit: 1): '))
    except ValueError:
        fitting = 0

    x1 = timeAxis[0]
    x1_scaled = timeAxis[1]
    idx1 = np.abs(x1_scaled - 1).argmin()
    x1_label = r'$t$/$\tau$'
    x2 = timeAxis[2]
    x2_scaled = timeAxis[3]
    idx2 = np.abs(x2_scaled - 0).argmin()
    x2_label = r'log[$t$/$\tau$]'

    if select == 0:
        y1 = creepComp(infMod, k, retardTime, x1)
        y1 *= 10**6             # rescale to MPa
        y1_label = r'$J$($t$) /MPa$^{{-1}}$'
        label1 = r'Creep compliance (linear)'
        creepComp = creepComp(infMod, k, retardTime, x2)
        y2 = np.log10(creepComp)                 # creepCompの出力をnp.ndarrayに変更したのでこのように簡単に書ける
        y2_label = r'log[ $J$($t$) /Pa$^{{-1}}$]'
        y_lim = [np.min(y2)-0.2, np.max(y2)+0.2]
        label2 = r'Creep compliance (log)'
        if fitting == 1:
            pass
        if fitting == 0:
            fitTimes = fitTimes()
            y_fit, param = curveFit(x2_scaled, y2, fitTimes)        
            fit_result = r'$J(t) \propto (t/\tau)^{{{0:.2f}}}$'.format(param[0])
            label_fit = r'Fitted creep compliance'
        legend_loc='upper left'
        savefile = './png/SLS1_creep_compliance.png'

    if select == 1:
        y1 = creepFunc(k, retardTime, x1)
        y1_label = r'$\psi$($t$)'
        label1 = 'Creep function (linear)'
        creepFunc = creepFunc(k, retardTime, x2)
        y2 = np.log10(creepFunc)                 # creepFuncの出力をnp.ndarrayに変更したのでこのように簡単に書ける
        y2_label = r'log[$\psi$($t$)]'
        y_lim = [np.min(y2)-1.0, np.max(y2)+1.0]
        label2 = 'Creep function (log)'
        if fitting == 1:
            pass
        if fitting == 0:
            fitTimes = fitTimes()
            y_fit, param = curveFit(x2_scaled, y2, fitTimes)        
            fit_result = r'$\psi(t) \propto (t/\tau)^{{{0:.2f}}}$'.format(param[0])
            label_fit = r'Fitted creep function'
        legend_loc='upper left'
        savefile = './png/SLS1_creep_func.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,10), tight_layout=True)
    ax1 = fig.add_subplot(211)
    ax1.set_title('SLS I model '+param_text)
    ax1.set_xlabel(x1_label)
    ax1.set_ylabel(y1_label)
    ax1.set_xlim(x1_scaled[0], x1_scaled[-1])
    ax1.set_ylim(-np.max(y1)*0.05, np.max(y1)*1.1)
    ax1.grid()
    ax1.plot(x1_scaled, y1, c='r', label=label1)
    ax1.legend(loc=legend_loc)
    ax1.vlines(1, 0, y1[idx1], color='k', ls='dashed')
    ax1.set_axisbelow(True)
    ax2 = fig.add_subplot(212)
    ax2.set_xlabel(x2_label)
    ax2.set_ylabel(y2_label)
    ax2.set_xlim(x2_scaled[0], x2_scaled[-1])
    ax2.set_ylim(y_lim[0], y_lim[1])
    ax2.grid()
    ax2.plot(x2_scaled, y2, c='b', label=label2)
    ax2.legend(loc=legend_loc)
    ax2.vlines(0, 0, y2[idx2], color='k', ls='dashed')
    ax2.set_axisbelow(True)
    if fitting == 0:
        ax2.plot(x2_scaled, y_fit, c='b', ls=':', label=label_fit)
        fig.text(0.7, 0.30, fit_result)
    
    fig.savefig(savefile, dpi=300)

    plt.show()