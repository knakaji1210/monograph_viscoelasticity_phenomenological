import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_complexComp(E, T, f):
    numer = E**(-1)
    denom = 1 + T*f*(2j/2)
    comComp = numer/denom
    strComp = comComp.real
    losComp = -comComp.imag
    return strComp, losComp

def complexComp(E, T, freq):
    strComp = []
    losComp = []
    for f in freq:
        s, l = calc_complexComp(E, T, f)
        strComp.append(s)
        losComp.append(l)
    complexComp = [strComp, losComp]
    return complexComp

def reqFreqs():
    try:
        minFreq = int(input('Enter minimum frequency in log scale (default = -1): '))
    except ValueError:
        minFreq = -1
    try:
        maxFreq = int(input('Enter maximum frequency in log scale (default = 5): '))
    except ValueError:
        maxFreq = 5
    try:
        minFitFreq1 = float(input('Enter minimum frequency for fitting (storage) (default = -1): '))
    except ValueError:
        minFitFreq1 = -1
    try:
        maxFitFreq1 = float(input('Enter maximum frequency for fitting (storage) (default = 5): '))
    except ValueError:
         maxFitFreq1 = 5
    try:
        minFitFreq2 = float(input('Enter minimum frequency for fitting (loss) (default = -1): '))
    except ValueError:
        minFitFreq2 = -1
    try:
        maxFitFreq2 = float(input('Enter maximum frequency for fitting (loss) (default = 5): '))
    except ValueError:
         maxFitFreq2 = 5
    intFreq = maxFreq - minFreq
    freqs = [minFreq, maxFreq, intFreq, minFitFreq1, maxFitFreq1, minFitFreq2, maxFitFreq2]
    return freqs

def reqParams():
    try:
        modulus = float(input('Enter modulus value (MPa) (default = 1 MPa): '))*10**6
    except ValueError:
        modulus = 10**6
    try:
        viscosity = float(input('Enter viscosity value (kPa s) (default = 10 kPa s): '))*10**3
    except ValueError:
        viscosity = 10**4
    relaxTime = viscosity/modulus
    return modulus, relaxTime

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
    E, T = reqParams()
    freqs = reqFreqs()
    freq = np.logspace(freqs[0], freqs[1], freqs[2]*4+1)
    strComp = complexComp(E, T, freq)[0]
    losComp = complexComp(E, T, freq)[1]
    log_freq = [np.log10(f) for f in freq]
    log_strComp = [np.log10(s) for s in strComp]
    log_losComp = [np.log10(l) for l in losComp]

    param_text = '(E = {0:.1f} MPa, T = {1:.1f} ms)'.format(E/10**6, T*10**3)

    # fitting
    minFit1 = fitRegion(freq, freqs[3], freqs[4])[0]
    maxFit1 = fitRegion(freq, freqs[3], freqs[4])[1]
    minFit2 = fitRegion(freq, freqs[5], freqs[5])[0]
    maxFit2 = fitRegion(freq, freqs[6], freqs[6])[1]
    paramStr,_ = curve_fit(loglogFit, log_freq[minFit1:maxFit1], log_strComp[minFit1:maxFit1])
    paramLos,_ = curve_fit(loglogFit, log_freq[minFit2:maxFit2], log_losComp[minFit2:maxFit2])
    fit_strComp = fittedArray(log_freq, paramStr)
    fit_losComp = fittedArray(log_freq, paramLos)
    fit_result1 = 'slope of storage compliance: {0:.2f}'.format(paramStr[0])
    fit_result2 = 'slope of loss compliance: {0:.2f}'.format(paramLos[0])

    # drawing graphs
    fig = plt.figure(figsize=(8,5))
    fig.suptitle('Kelvin-Voigt model '+param_text)
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('log(frequency /s)')
    ax1.set_ylabel('log(compliance /Pa^(-1))')
    ax1.set_xlim(freqs[0],freqs[1])
    ax1.set_ylim(np.min(log_losComp)*1.2, np.max(log_losComp)*0.8)
    ax1.grid()
    ax1.scatter(log_freq, log_strComp, c='r', label='Storage Compliance')
    ax1.scatter(log_freq, log_losComp, c='b', label='Loss Compliance')
    ax1.plot(log_freq, fit_strComp, c='r', ls=':', label='fitted Storage Compliance')
    ax1.plot(log_freq, fit_losComp, c='b', ls=':', label='fitted Loss Compliance')
    ax1.legend(loc='lower left')

    fig.text(0.6, 0.75, fit_result1)
    fig.text(0.6, 0.70, fit_result2)

    fig.savefig('./png/voigt_model_complex.png')

    plt.show()