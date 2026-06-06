import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_relaxFunc(E, k, T, t):
    relaxFunc = (1 - 1/k)*(1 - np.exp(-t/T))
    return relaxFunc

def relaxFunc(E, k, T, tim):
    relaxFunc = []
    for t in tim:
        r = calc_relaxFunc(E, k, T, t)
        relaxFunc.append(r)
    return relaxFunc

def reqTimes():
    try:
        minTime = int(input('Enter minimum time in log scale (default = -6): '))
    except ValueError:
        minTime = -6
    try:
        maxTime = int(input('Enter maximum time in log scale (default = 0): '))
    except ValueError:
        maxTime = 0
    try:
        minFitTime = float(input('Enter minimum time for fitting (default = -4): '))
    except ValueError:
        minFitTime = -4
    try:
        maxFitTime = float(input('Enter maximum time for fitting (default = -3): '))
    except ValueError:
         maxFitTime = -3
    intTime = maxTime - minTime
    Times = [minTime, maxTime, intTime, minFitTime, maxFitTime]
    return Times

def reqParams():
    try:
        insMod = float(input('Enter instantaneous modulus value (GPa) (default = 1 GPa): '))*10**9
    except ValueError:
        insMod = 10**9
    try:
        infMod = float(input('Enter equilibrium modulus value (MPa) (default = 0.1 MPa): '))*10**5
    except ValueError:
        infMod = 10**5
    try:
        viscosity = float(input('Enter viscosity value (kPa s) (default = 100 kPa s): '))*10**5
    except ValueError:
        viscosity = 10**5
    modulus = insMod - infMod
    relaxTime = viscosity/modulus
    k = insMod/infMod
    return insMod, k, relaxTime

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
    # calculating Relaxation modulus
    insMod, k, relaxTime = reqParams()
    Times = reqTimes()
    tim = np.logspace(Times[0], Times[1], Times[2]*6+1)
    relaxFunc = relaxFunc(insMod, k, relaxTime, tim)
    log_time = [np.log10(t) for t in tim]
    log_relaxFunc = [np.log10(r) for r in relaxFunc]

    param_text = '(Ei = {0:.1f} MPa, Einf = {1:.1f} MPa, T = {2:.1f} ms)'.format(insMod/10**6, insMod/(k*10**6), relaxTime*10**3)

    # fitting
    minFit = fitRegion(tim, Times[3], Times[4])[0]
    maxFit = fitRegion(tim, Times[3], Times[4])[1]
    paramComp,_ = curve_fit(loglogFit, log_time[minFit:maxFit], log_relaxFunc[minFit:maxFit])
    fit_relaxFunc = fittedArray(log_time, paramComp)
    fit_result = 'slope of relaxation function: {0:.2f}'.format(paramComp[0])

    # drawing graphs
    fig = plt.figure(figsize=(8,5))
    fig.suptitle('SLS model II '+param_text)
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('log(Time /s)')
    ax1.set_ylabel('log(Relaxation function)')
    ax1.set_xlim(Times[0],Times[1])
    ax1.set_ylim(np.min(log_relaxFunc)-1, np.max(log_relaxFunc)+1)
    ax1.grid()
    ax1.scatter(log_time, log_relaxFunc, c='r', label='Relaxation function')
    ax1.plot(log_time, fit_relaxFunc, c='b', ls=':', label='fitted Relaxation funcion')
    ax1.legend(loc='upper left')

    fig.text(0.4, 0.5, fit_result)

    fig.savefig('./png/sls_model_relax_func.png')

    plt.show()