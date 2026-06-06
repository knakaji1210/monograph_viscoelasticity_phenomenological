import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_complexMod(k0, k, c, f):
# complexを使うバージョン
#    numer = complex(0, k*c*f)
#    denom = complex(k, c*f)
# (2j/2)を使うバージョン
    numer = k*c*f*(2j/2)
    denom = k + c*f*(2j/2)
    comMod = k0 + numer/denom
    strMod = comMod.real
    losMod = comMod.imag
    return strMod, losMod

def complexMod(k0, k, c, freq):
    strMod = []
    losMod = []
    for f in freq:
        s, l = calc_complexMod(k0, k, c, f)
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
    k_list = []
    c_list = []
    try:
        k0 = float(input('Enter k0 (default = 1): '))
    except ValueError:
        k0 = 1
    try:
        num = int(input('Enter the number of Maxwell components (default = 1): '))
    except ValueError:
        num = 1   
    for i in range(num):
        try:
            k = float(input('Enter k (default = 10): '))
        except ValueError:
            k = 10
        k_list.append(k)
        try:
            c = float(input('Enter c (default = 100): '))
        except ValueError:
            c = 100
        c_list.append(c)
    params = [num, k0, k_list, c_list]
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
    strMod = [params[1] for s in range(len(freq))]
    losMod = [0 for l in range(len(freq))]
    for i in range(params[0]):
        # k0は系で１つのバージョン
        strMod = [s1 + s2 for (s1, s2) in zip(strMod, complexMod(0, params[2][i], params[3][i], freq)[0])]
        losMod = [l1 + l2 for (l1, l2) in zip(losMod, complexMod(0, params[2][i], params[3][i], freq)[1])]
    strMod = [s + params[1] for s in strMod]
    log_freq = [np.log10(f) for f in freq]
    log_strMod = [np.log10(s) for s in strMod]
    log_losMod = [np.log10(l) for l in losMod]
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
    plt.title('Standard Linear Viscoelastic Model (Multiple)')
    plt.xlabel('Frequency')
    plt.ylabel('Modulus')
    ax1 = plt.scatter(log_freq,log_strMod, c='r', label='Storage Modulus')
    ax2 = plt.scatter(log_freq,log_losMod, c='b', label='Loss Modulus')
    ax3 = plt.plot(log_freq, fit_strMod, c='r', ls=':', label='fitted Storage Modulus')
    ax4 = plt.plot(log_freq, fit_losMod, c='b', ls=':', label='fitted Loss Modulus')
    plt.legend(loc='upper left')
    plt.xlim(freqs[0],freqs[1])
    plt.ylim(ymin-1,ymax+3)
    k0_text = 'k0 = {0:.1f}'.format(params[1])
    plt.text(freqs[0]+0.2, (ymin+ymax)/2 + 2.5, k0_text)   
    for i in range(params[0]):
        param_text = 'k{0} = {1:.1f}, c{0} = {2:.1f}'.format(i+1, params[2][i], params[3][i])
        plt.text(freqs[0]+0.2, (ymin+ymax)/2 + 2.0 - 0.5*i, param_text)
    plt.text((freqs[0]+freqs[1])/2, (ymin+ymax)/2 - 2, fit_result1)
    plt.text((freqs[0]+freqs[1])/2, (ymin+ymax)/2 - 2.5, fit_result2)

    fig.savefig('./png/slv_model_complex.png')

    plt.show()