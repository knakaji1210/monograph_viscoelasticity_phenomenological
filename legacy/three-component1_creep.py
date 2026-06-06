import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_creepFunc(E, k, T, t):
    creepFunc = (1 - (1 - 1/k)*np.exp(-t/T))/E
    return creepFunc

def creepFunc(E, k, T, tim):
    creepFunc = []
    for t in tim:
        c = calc_creepFunc(E, k, T, t)
        creepFunc.append(c)
    return creepFunc

def reqTimes():
    try:
        minTime = int(input('Enter minimum time in log scale (default = -8): '))
    except ValueError:
        minTime = -8
    try:
        maxTime = int(input('Enter maximum time in log scale (default = 2): '))
    except ValueError:
        maxTime = 2
    try:
        minFitTime = float(input('Enter minimum time for fitting (default = -4): '))
    except ValueError:
        minFitTime = -4
    try:
        maxFitTime = float(input('Enter maximum time for fitting (default = -2): '))
    except ValueError:
         maxFitTime = -2
    intTime = maxTime - minTime
    Times = [minTime, maxTime, intTime, minFitTime, maxFitTime]
    return Times

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
        viscosity = float(input('Enter viscosity value (kPa s) (default = 10 kPa s): '))*10**3
    except ValueError:
        viscosity = 10**4
    modulus = insMod*infMod/(insMod - infMod)
    relaxTime = viscosity/modulus
    k = insMod/infMod
    return infMod, k, relaxTime

def getNearestIdx(list, num):
    idx = np.abs(np.asarray(list) - num).argmin()
    return idx

def fitRegion(tim, minNum, maxNum):
    minFit = getNearestIdx(tim, 10**minNum)
    maxFit = getNearestIdx(tim, 10**maxNum)
    fitRegion = [minFit, maxFit]
    return fitRegion

def loglogFit(x, a, b):
    return  a*x + b

def fittedArray(x_array, param):
    fitted_array = [loglogFit(num, param[0], param[1]) for num in x_array]
    return fitted_array

if __name__=='__main__':
    # calculating Creep compliance
    infMod, k, relaxTime = reqParams()
    Times = reqTimes()
    tim = np.logspace(Times[0], Times[1], Times[2]*4+1)
    creepFunc = creepFunc(infMod, k, relaxTime, tim)
    log_time = [np.log10(t) for t in tim]
    log_creepFunc = [np.log10(c) for c in creepFunc]

    param_text = '(Ei = {0:.1f} MPa, Einf = {1:.1f} MPa, T = {2:.1f} ms)'.format(infMod*k/10**6, infMod/10**6, relaxTime*10**3)

    # fitting
    minFit = fitRegion(tim, Times[3], Times[4])[0]
    maxFit = fitRegion(tim, Times[3], Times[4])[1]
    paramComp,_ = curve_fit(loglogFit, log_time[minFit:maxFit], log_creepFunc[minFit:maxFit])
    fit_creepFunc = fittedArray(log_time, paramComp)
    fit_result = 'slope of creep compliance: {0:.2f}'.format(paramComp[0])

    # drawing graphs
    fig = plt.figure(figsize=(8,5))
    fig.suptitle('SLS model I '+param_text)
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('log(Time /s)')
    ax1.set_ylabel('log(Creep compliance /Pa^(-1))')
    ax1.set_xlim(Times[0],Times[1])
    ax1.set_ylim(np.min(log_creepFunc)-1, np.max(log_creepFunc)+1)
    ax1.grid()
    ax1.scatter(log_time, log_creepFunc, c='r', label='Creep compliance')
    ax1.plot(log_time, fit_creepFunc, c='b', ls=':', label='fitted Creep compliance')
    ax1.legend(loc='upper left')

    fig.text(0.6, 0.3, fit_result)

    fig.savefig('./png/sls_model_creep.png')

    plt.show()