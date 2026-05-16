# relaxation modulus of SLS model 2

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_relaxMod(E, k, T, t):
    relaxMod = E*(1/k + (1 - 1/k)*np.exp(-t/T))
    # E must be insMod
    return relaxMod

def relaxMod(E, k, T, tim):
    relaxMod = []
    for t in tim:
        c = calc_relaxMod(E, k, T, t)
        relaxMod.append(c)
    return relaxMod

def calc_relaxFunc(T, t):
    relaxFunc = (1 - 1/k)*(1 - np.exp(-t/T))
    return relaxFunc

def relaxFunc(T, tim):
    relaxFunc = []
    for t in tim:
        r = calc_relaxFunc(T, t)
        relaxFunc.append(r)
    return relaxFunc

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
    relaxTime = viscosity/modulus
    k = insMod/infMod
    return insMod, infMod, k, relaxTime

def timeAxis(relaxTime):
    log_relaxT = np.log10(relaxTime)
    Time = np.logspace(log_relaxT-2, log_relaxT+3, 31)
    scaledTime = [np.log10(t/relaxTime) for t in Time]
    timeAxis = [Time, scaledTime]
    return timeAxis

def fitTimes():
    try:
        minTime = float(input('Enter minimum time for fitting in log scale (default = 0): '))
    except ValueError:
        minTime = 0
    try:
        maxTime = float(input('Enter maximum time for fitting in log scale (default = 0.5): '))
    except ValueError:
         maxTime = 0.5
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
    # calculating relaxation modulus and relaxation funcion
    insMod, infMod, k, relaxTime = reqParams()
    param_text = '(Ei = {0:.1f} MPa, Einf = {1:.1f} MPa, tau = {2:.1f} ms)'.format(insMod/10**6, infMod/10**6, relaxTime*10**3)
    timeAxis = timeAxis(relaxTime)
    try:
        select = int(input('Selection (relaxation modulus: 0, relaxation function: 1): '))
    except ValueError:
        select = 0
    
    if select == 0:
        x1 = timeAxis[0]
        x1_scaled = timeAxis[1]
        x1_label = 'log[t/tau]'
        relaxMod = relaxMod(insMod, k, relaxTime, x1)  # E = insMod
        y1 = [np.log10(r) for r in relaxMod]
        y1_label = 'log[E(t) /Pa]'
        label1 = 'Relaxation modulus (log)'
        fitTimes = fitTimes()
        y2, param = curveFit(x1_scaled, y1, fitTimes)
        label2 = 'fitted relaxation modulus'
        fit_result = 'E ∝ (t/tau)^({0:.2f})'.format(param[0])
        legend_loc='upper right'
        savefile = './png/SLS2_relaxation_modulus.png'

    if select == 1:
        x1 = timeAxis[0]
        x1_scaled = timeAxis[1]
        x1_label = 'log[t/tau]'
        relaxFunc = relaxFunc(relaxTime, x1)
        y1 = [np.log10(r) for r in relaxFunc]
        y1_label = 'log[phi(t)]'
        label1 = 'Relaxation function (log)'
        fitTimes = fitTimes()
        y2, param = curveFit(x1_scaled, y1, fitTimes)
        label2 = 'fitted relaxation function'
        fit_result = 'E ∝ (t/tau)^({0:.2f})'.format(param[0])
        legend_loc='upper left'
        savefile = './png/SLS2_relaxation_func.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title('SLS model II '+param_text)
    ax.set_xlabel(x1_label)
    ax.set_ylabel(y1_label)
    ax.set_xlim(-2.2,3.2)
    ax.set_ylim(np.min(y1)-0.5, np.max(y1)+0.5)
    ax.grid()
    ax.scatter(x1_scaled, y1, c='r', label=label1)
    ax.plot(x1_scaled, y2, c='b', ls=':', label=label2)
    ax.legend(loc=legend_loc)
    ax.set_axisbelow(True)
    fig.text(0.55, 0.5, fit_result)

    fig.savefig(savefile)

    plt.show()