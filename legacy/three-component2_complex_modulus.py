import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_complexMod(E, k, T, f):
    numer = 1/k + T*f*(2j/2)
    denom = 1 + T*f*(2j/2)
    comMod = E*numer/denom
    strMod = comMod.real
    losMod = comMod.imag
    return strMod, losMod

def complexMod(E, k, T, freq):
    strMod = []
    losMod = []
    for f in freq:
        s, l = calc_complexMod(E, k, T, f)
        strMod.append(s)
        losMod.append(l)
    dynamicMod = [strMod, losMod]
    return dynamicMod

def reqFreqs():
    try:
        minFreq = int(input('Enter minimum frequency in log scale (default = -2): '))
    except ValueError:
        minFreq = -2
    try:
        maxFreq = int(input('Enter maximum frequency in log scale (default = 6): '))
    except ValueError:
        maxFreq = 6
    try:
        minFitFreq1 = float(input('Enter minimum frequency for fitting (storage) (default = 1.2): '))
    except ValueError:
        minFitFreq1 = 1.2
    try:
        maxFitFreq1 = float(input('Enter maximum frequency for fitting (storage) (default = 2.7): '))
    except ValueError:
         maxFitFreq1 = 2.7
    try:
        minFitFreq2 = float(input('Enter minimum frequency for fitting (loss) (default = -2): '))
    except ValueError:
        minFitFreq2 = -2
    try:
        maxFitFreq2 = float(input('Enter maximum frequency for fitting (loss) (default = 2): '))
    except ValueError:
         maxFitFreq2 = 2
    intFreq = maxFreq - minFreq
    freqs = [minFreq, maxFreq, intFreq, minFitFreq1, maxFitFreq1, minFitFreq2, maxFitFreq2]
    return freqs

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
        viscosity = float(input('Enter viscosity value (kPa s) (default = 1 MPa s): '))*10**6
    except ValueError:
        viscosity = 10**6
    modulus = insMod - infMod
    relaxTime = viscosity/modulus
    k = insMod/infMod
    return insMod, k, relaxTime

def getNearestIdx(list, num):
    idx = np.abs(np.asarray(list) - num).argmin()
    return idx

def fitRegion(freq, minNum, maxNum):
    minFit = getNearestIdx(freq, 10**minNum)
    maxFit = getNearestIdx(freq, 10**maxNum)
    fitRegion = [minFit, maxFit]
    return fitRegion

def loglogFit(x, a, b):
    return  a*x + b

def fittedArray(x_array, param):
    fitted_array = [loglogFit(num, param[0], param[1]) for num in x_array]
    return fitted_array

if __name__=='__main__':
    # calculating dynamic Modulus
    insMod, k, relaxTime = reqParams()
    freqs = reqFreqs()
    freq = np.logspace(freqs[0], freqs[1], freqs[2]*4+1)
    strMod = complexMod(insMod, k, relaxTime, freq)[0]
    losMod = complexMod(insMod, k, relaxTime, freq)[1]
    log_freq = [np.log10(f) for f in freq]
    log_strMod = [np.log10(s) for s in strMod]
    log_losMod = [np.log10(l) for l in losMod]

    param_text = '(Ei = {0:.1f} GPa, Einf = {1:.1f} MPa, T = {2:.1f} ms)'.format(insMod/10**9, insMod/(k*10**6), relaxTime*10**3)

    # fitting
    minFit1 = fitRegion(freq, freqs[3], freqs[4])[0]
    maxFit1 = fitRegion(freq, freqs[3], freqs[4])[1]
    minFit2 = fitRegion(freq, freqs[5], freqs[5])[0]
    maxFit2 = fitRegion(freq, freqs[6], freqs[6])[1]
    paramStr,_ = curve_fit(loglogFit, log_freq[minFit1:maxFit1], log_strMod[minFit1:maxFit1])
    paramLos,_ = curve_fit(loglogFit, log_freq[minFit2:maxFit2], log_losMod[minFit2:maxFit2])
    fit_strMod = fittedArray(log_freq, paramStr)
    fit_losMod = fittedArray(log_freq, paramLos)
    fit_result1 = 'slope of storage modulus: {0:.2f}'.format(paramStr[0])
    fit_result2 = 'slope of loss modulus: {0:.2f}'.format(paramLos[0])

    # drawing graphs
    fig = plt.figure(figsize=(8,5))
    fig.suptitle('SLS model II '+param_text)
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('log(frequency /s)')
    ax1.set_ylabel('log(modulus /Pa)')
    ax1.set_xlim(freqs[0],freqs[1])
    ax1.set_ylim(np.min(log_losMod)-1, np.max(log_losMod)+1)
    ax1.grid()
    ax1.scatter(log_freq, log_strMod, c='r', label='Storage Modulus')
    ax1.scatter(log_freq, log_losMod, c='b', label='Loss Modulus')
    ax1.plot(log_freq, fit_strMod, c='r', ls=':', label='fitted Storage Modulus')
    ax1.plot(log_freq, fit_losMod, c='b', ls=':', label='fitted Loss Modulus')
    ax1.legend(loc='upper left')

    fig.text(0.6, 0.30, fit_result1)
    fig.text(0.6, 0.25, fit_result2)

    fig.savefig('./png/sls_model_complex.png')

    plt.show()