import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_complexMod(a, b, c, beta, f):
    numer = a + b * (2j/2)**beta * f**beta
    denom = 1 + c * (2j/2)**beta * f**beta
    comMod = numer/denom
    strMod = comMod.real
    losMod = comMod.imag
    return strMod, losMod

def complexMod(a, b, c, beta, freq):
    strMod = []
    losMod = []
    for f in freq:
        s, l = calc_complexMod(a, b, c, beta, f)
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
    a_list = []
    b_list = []
    c_list = []
    beta_list = []
    try:
        num = int(input('Enter the number of fractional components (default = 1): '))
    except ValueError:
        num = 1   
    for i in range(num):
        try:
            a = float(input('Enter a (default = 200): '))
        except ValueError:
            a = 200
        a_list.append(a)
        try:
            b = float(input('Enter b (default = 500): '))
        except ValueError:
            b = 500
        b_list.append(b)
        try:
            c = float(input('Enter c (default = 0.01): '))
        except ValueError:
            c = 0.01
        c_list.append(c)
        try:
            beta = float(input('Enter beta (default = 1): '))
        except ValueError:
            beta = 1
        beta_list.append(beta)
    params = [num, a_list, b_list, c_list, beta_list]
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
    strMod = [0 for s in range(len(freq))]
    losMod = [0 for l in range(len(freq))]
    for i in range(params[0]):
        strMod = [s1 + s2 for (s1, s2) in zip(strMod, complexMod(params[1][i], params[2][i], params[3][i], params[4][i], freq)[0])]
        losMod = [l1 + l2 for (l1, l2) in zip(losMod, complexMod(params[1][i], params[2][i], params[3][i], params[4][i], freq)[1])]
    log_freq = [np.log10(f) for f in freq]
    log_strMod = [np.log10(s) for s in strMod]
    log_losMod = [np.log10(l) for l in losMod]
    ymin = np.min([log_strMod[0],log_losMod[0]])
    ymax = log_strMod[-1]
    print(params[4])

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
    plt.title('Fractional Calculus Model (Multiple)')
    plt.xlabel('Frequency')
    plt.ylabel('Modulus')
    ax1 = plt.scatter(log_freq,log_strMod, c='r', label='Storage Modulus')
    ax2 = plt.scatter(log_freq,log_losMod, c='b', label='Loss Modulus')
    ax3 = plt.plot(log_freq, fit_strMod, c='r', ls=':', label='fitted Storage Modulus')
    ax4 = plt.plot(log_freq, fit_losMod, c='b', ls=':', label='fitted Loss Modulus')
    plt.legend(loc='upper left')
    plt.xlim(freqs[0],freqs[1])
    plt.ylim(ymin-1,ymax+3)
    for i in range(params[0]):
        param_text1 = 'a{0} = {1:.1f}, b{0} = {2:.1f}'.format(i+1, params[1][i], params[2][i])
        param_text2 = 'c{0} = {1:.2f}, beta{0} = {2:.2f}'.format(i+1, params[3][i], params[4][i])
        plt.text(freqs[0]+0.2, (ymin+ymax)/2 + 2.0 - i, param_text1)
        plt.text(freqs[0]+0.2, (ymin+ymax)/2 + 1.5 - i, param_text2)
    plt.text((freqs[0]+freqs[1])/2, (ymin+ymax)/2 - 2, fit_result1)
    plt.text((freqs[0]+freqs[1])/2, (ymin+ymax)/2 - 2.5, fit_result2)

    fig.savefig('./png/frac_model.png')

    plt.show()