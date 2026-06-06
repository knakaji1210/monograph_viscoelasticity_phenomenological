import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_dynamicMod(k0, k, c, f):
    strMod = k0 + k * c**2 * f**2 / (k**2 + c**2 * f**2)
    losMod = k**2 * c * f / (k**2 + c**2 * f**2)
    return strMod, losMod

def dynamicMod(k0, k, c, freq):
    strMod = []
    losMod = []
    for f in freq:
        s, l = calc_dynamicMod(k0, k, c, f)
        strMod.append(s)
        losMod.append(l)
    dynamicMod = [strMod, losMod]
    return dynamicMod

def reqFreqs():
    try:
        minFreq = int(input('Enter minimum frequency in log scale (default = -3): '))
    except ValueError:
        minFreq = -3
    try:
        maxFreq = int(input('Enter maximum frequency in log scale (default = 6): '))
    except ValueError:
        maxFreq = 6
    try:
        minFitFreq1 = float(input('Enter minimum frequency for fitting (storage) (default = -3): '))
    except ValueError:
        minFitFreq1 = -3
    try:
        maxFitFreq1 = float(input('Enter maximum frequency for fitting (storage) (default = 6): '))
    except ValueError:
         maxFitFreq1 = 6
    try:
        minFitFreq2 = float(input('Enter minimum frequency for fitting (loss) (default = -3): '))
    except ValueError:
        minFitFreq2 = -3
    try:
        maxFitFreq2 = float(input('Enter maximum frequency for fitting (loss) (default = 6): '))
    except ValueError:
         maxFitFreq2 = 6
    intFreq = maxFreq - minFreq
    freqs = [minFreq, maxFreq, intFreq, minFitFreq1, maxFitFreq1, minFitFreq2, maxFitFreq2]
    return freqs

def reqParams():
    try:
        k0 = float(input('Enter k0 (default = 1): '))
    except ValueError:
        k0 = 1
    try:
        k = float(input('Enter k (default = 10): '))
    except ValueError:
        k = 10
    try:
        c = float(input('Enter c (default = 100): '))
    except ValueError:
        c = 100
    params = [k0, k, c]
    return params

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
    freqs = reqFreqs()
    freq = np.logspace(freqs[0], freqs[1], freqs[2]*4+1)
    params = reqParams()
    strMod = dynamicMod(params[0], params[1], params[2], freq)[0]
    losMod = dynamicMod(params[0], params[1], params[2], freq)[1]
    log_freq = [np.log10(f) for f in freq]
    log_strMod = [np.log10(s) for s in strMod]
    log_losMod = [np.log10(l) for l in losMod]
    param_text1 = 'k0 = {0:.1f}'.format(params[0])
    param_text2 = 'k = {0:.1f}, c = {1:.1f}'.format(params[1], params[2])
    ymin = np.min([log_strMod[0],log_losMod[0]])
    ymax = log_strMod[-1]

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
    fig = plt.figure()
    plt.title('Standard Linear Viscoelastic Model')
    plt.xlabel('Frequency')
    plt.ylabel('Modulus')
    ax1 = plt.scatter(log_freq,log_strMod, c='r', label='Storage Modulus')
    ax2 = plt.scatter(log_freq,log_losMod, c='b', label='Loss Modulus')
    ax3 = plt.plot(log_freq, fit_strMod, c='r', ls=':', label='fitted Storage Modulus')
    ax4 = plt.plot(log_freq, fit_losMod, c='b', ls=':', label='fitted Loss Modulus')
    plt.legend(loc='upper left')
    plt.xlim(freqs[0],freqs[1])
    plt.ylim(ymin-1,ymax+3)
    plt.text(freqs[0]+0.5, (ymin+ymax)/2 + 2.5, param_text1)
    plt.text(freqs[0]+0.5, (ymin+ymax)/2 + 2.0, param_text2)
    plt.text((freqs[0]+freqs[1])/2, (ymin+ymax)/2 - 2, fit_result1)
    plt.text((freqs[0]+freqs[1])/2, (ymin+ymax)/2 - 2.5, fit_result2)
    plt.show()