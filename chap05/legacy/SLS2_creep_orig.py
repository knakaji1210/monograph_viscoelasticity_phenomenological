# creep compliance of SLS model 2

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_creepComp(E, k, T, t):
    creepComp = (1 - (1 - 1/k)*np.exp(-t/(k*T)))/E
    # E must be infMod
    return creepComp

def creepComp(E, k, T, tim):
    creepComp = []
    for t in tim:
        c = calc_creepComp(E, k, T, t)
        creepComp.append(c)
    return creepComp

def calc_creepFunc(T, t):
    creepFunc = (k - 1)*(1 - np.exp(-t/T))
    return creepFunc

def creepFunc(T, tim):
    creepFunc = []
    for t in tim:
        r = calc_creepFunc(T, t)
        creepFunc.append(r)
    return creepFunc

def reqParams():
    try:
        insMod = float(input('Enter instantaneous modulus value (MPa) (default = 10 MPa): '))*10**6
    except ValueError:
        insMod = 10**7
    try:
        infMod = float(input('Enter equilibrium modulus value (MPa) (default = 1 MPa): '))*10**6
    except ValueError:
        infMod = 10**6
    try:
        viscosity = float(input('Enter viscosity value (kPa s) (default = 900 kPa s): '))*10**3
    except ValueError:
        viscosity = 9 * 10**5
    modulus = insMod - infMod
    retardTime = viscosity/modulus
    k = insMod/infMod
    return insMod, infMod, k, retardTime

def timeAxis(retardTime):
    log_relaxT = np.log10(retardTime)
    Time = np.logspace(log_relaxT-3, log_relaxT+3, 31)
    scaledTime = [np.log10(t/retardTime) for t in Time]
    timeAxis = [Time, scaledTime]
    return timeAxis

def fitTimes():
    try:
        minTime = float(input('Enter minimum time for fitting in log scale (default = 0.2): '))
    except ValueError:
        minTime = 0.2
    try:
        maxTime = float(input('Enter maximum time for fitting in log scale (default = 0.8): '))
    except ValueError:
         maxTime = 0.8
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
    param_text = '(Ei = {0:.1f} MPa, Einf = {1:.1f} MPa, tau = {2:.1f} ms)'.format(insMod/10**6, infMod/10**6, retardTime*10**3)
    timeAxis = timeAxis(retardTime)
    try:
        select = int(input('Selection (creep compliance: 0, creep function: 1): '))
    except ValueError:
        select = 0
    
    if select == 0:
        x1 = timeAxis[0]
        x1_scaled = timeAxis[1]
        x1_label = 'log[t/tau]'
        creepComp = creepComp(infMod, k, retardTime, x1)  # E = infMod
        y1 = [np.log10(r) for r in creepComp]
        y1_label = 'log[J(t) /Pa^(-1)]'
        label1 = 'Creep compliance (log)'
        fitTimes = fitTimes()
        y2, param = curveFit(x1_scaled, y1, fitTimes)
        label2 = 'fitted creep compliance'
        fit_result = 'J ∝ (t/tau)^({0:.2f})'.format(param[0])
        legend_loc='upper left'
        savefile = './png/SLS2_creep_compliance.png'

    if select == 1:
        x1 = timeAxis[0]
        x1_scaled = timeAxis[1]
        x1_label = 'log[t/tau]'
        creepFunc = creepFunc(retardTime, x1)
        y1 = [np.log10(r) for r in creepFunc]
        y1_label = 'log[psi(t)]'
        label1 = 'Creep function (log)'
        fitTimes = fitTimes()
        y2, param = curveFit(x1_scaled, y1, fitTimes)
        label2 = 'fitted creep function'
        fit_result = 'J ∝ (t/tau)^({0:.2f})'.format(param[0])
        legend_loc='upper left'
        savefile = './png/SLS2_creep_func.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title('SLS model II '+param_text)
    ax.set_xlabel(x1_label)
    ax.set_ylabel(y1_label)
    ax.set_xlim(-3.1,6.1)
    ax.set_ylim(np.min(y1)-0.5, np.max(y1)+0.5)
    ax.grid()
    ax.scatter(x1_scaled, y1, c='r', label=label1)
    ax.plot(x1_scaled, y2, c='b', ls=':', label=label2)
    ax.legend(loc=legend_loc)
    ax.set_axisbelow(True)
    fig.text(0.7, 0.30, fit_result)

    fig.savefig(savefile)

    plt.show()