import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_std_model(k0, k, c, f):
    strMod = k0 + k * c**2 * f**2 / (k**2 + c**2 * f**2)
    losMod = k**2 * c * f / (k**2 + c**2 * f**2)
    return strMod, losMod

def std_model(k0, k, c, freq):
    strMod = []
    losMod = []
    for f in freq:
        s, l = calc_std_model(k0, k, c, f)
        strMod.append(s)
        losMod.append(l)
    log_strMod = [np.log10(s) for s in strMod]
    log_losMod = [np.log10(l) for l in losMod]
    return log_strMod, log_losMod

def loglogFit(x, a, b):
    return  a*x + b

def fittedArray(x_array, param):
    fitted_array = [loglogFit(num, param[0], param[1]) for num in x_array]
    return fitted_array

if __name__=='__main__':
    minFreq = int(input('Enter minimum frequency in log scale: '))
    maxFreq = int(input('Enter maximum frequency in log scale: '))
    freqInterval = maxFreq - minFreq
    freq = np.logspace(minFreq, maxFreq, freqInterval*4+1)
    log_freq = [np.log10(f) for f in freq]
    log_strMod, log_losMod = std_model(1, 100, 100, freq)

    # fitting
    minFit = 0
    maxFit = int(len(freq)/2)
    param,_ = curve_fit(loglogFit, log_freq[minFit:maxFit], log_losMod[minFit:maxFit])
    fit_losMod = fittedArray(log_freq, param)
    fit_result = 'slope: {0}'.format(param[0])

    fig = plt.figure()
    plt.title('Standard Linear Viscoelastic Model')
    plt.xlabel('Frequency')
    plt.ylabel('Modulus')
    ax1 = plt.scatter(log_freq,log_strMod, c='r', label='Storage Modulus')
    ax2 = plt.scatter(log_freq,log_losMod, c='b', label='Loss Modulus')
    ax3 = plt.plot(log_freq, fit_losMod)
    plt.legend(loc='upper left')
    plt.xlim(minFreq,maxFreq)
    plt.ylim(-2,6)
    plt.text(-2.5, 4, fit_result)
    plt.show()