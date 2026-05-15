# complex compliance of SLS model 1

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_complexComp(E, k, T, f):
    numer = (k + T*f*(2j/2))/E
    # E must be insMod
    denom = 1 + T*f*(2j/2)
    comComp = numer/denom
    strComp = comComp.real
    losComp = -comComp.imag
    return strComp, losComp

def complexComp(E, k, T, freq):
    strComp = []
    losComp = []
    for f in freq:
        s, l = calc_complexComp(E, k, T, f)
        strComp.append(s)
        losComp.append(l)
    complexComp = [strComp, losComp]
    return complexComp

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
        viscosity = float(input('Enter viscosity value (kPa s) (default = 111.1 kPa s): '))*10**3
    except ValueError:
        viscosity = 10**6 / 9
    modulus = insMod*infMod/(insMod - infMod)
    retardTime = viscosity/modulus
    k = insMod/infMod
    return insMod, infMod, k, retardTime

def freqAxis(relaxTime):
    centerFreq = 1 / relaxTime
    freq = np.logspace(int(np.log10(centerFreq))-3, int(np.log10(centerFreq))+3, 43)
    scaledFreq = [f*relaxTime for f in freq]
    freqAxis = [freq, scaledFreq]
    return freqAxis

def fitFreqs():
    try:
        minFreq1 = float(input('Enter minimum frequency for fitting (storage) (default = 0): '))
    except ValueError:
        minFreq1 = 0
    try:
        maxFreq1 = float(input('Enter maximum frequency for fitting (storage) (default = 0.5): '))
    except ValueError:
         maxFreq1 = 0.5
    try:
        minFreq2 = float(input('Enter minimum frequency for fitting (loss) (default = -3): '))
    except ValueError:
        minFreq2 = -3
    try:
        maxFreq2 = float(input('Enter maximum frequency for fitting (loss) (default = -1): '))
    except ValueError:
         maxFreq2 = -1
    fitFreqs = [minFreq1, maxFreq1, minFreq2, maxFreq2]
    return fitFreqs

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

def curveFit(freq, fitFreqs):
    minFit1 = fitRegion(scaledFreq, fitFreqs[0], fitFreqs[1])[0]
    maxFit1 = fitRegion(scaledFreq, fitFreqs[0], fitFreqs[1])[1]
    minFit2 = fitRegion(scaledFreq, fitFreqs[2], fitFreqs[3])[0]
    maxFit2 = fitRegion(scaledFreq, fitFreqs[2], fitFreqs[3])[1]
    paramStr,_ = curve_fit(loglogFit, log_freq[minFit1:maxFit1], log_strComp[minFit1:maxFit1])
    paramLos,_ = curve_fit(loglogFit, log_freq[minFit2:maxFit2], log_losComp[minFit2:maxFit2])
    fit_strComp = fittedArray(log_freq, paramStr)
    fit_losComp = fittedArray(log_freq, paramLos)
    fit_result1 = "J' ∝ (freq*tau)^({0:.2f})".format(paramStr[0])
    fit_result2 = 'J" ∝ (freq*tau)^({0:.2f})'.format(paramLos[0])
    return fit_strComp, fit_losComp, fit_result1, fit_result2

if __name__=='__main__':
    # calcul1ating complex compliance and loss tangent
    insMod, infMod, k, relaxTime = reqParams()
    param_text = '(Ei = {0:.1f} MPa, Einf = {1:.1f} MPa, tau = {2:.1f} ms)'.format(insMod/10**6, infMod/10**6, relaxTime*10**3)
    freqAxis = freqAxis(relaxTime)
    fitting = -1
    try:
        select = int(input('Selection (complex compliance (linear): 0, complex compliance (log): 1, loss tangent: 2): '))
    except ValueError:
        select = 0

    freq = freqAxis[0]
    scaledFreq = freqAxis[1]
    freq_label = 'log[freq*tau]'
    strComp = complexComp(insMod, k, relaxTime, freq)[0]  # E = insMod
    losComp = complexComp(insMod, k, relaxTime, freq)[1]  # E = insMod
    losTan = [l/s for (s, l) in zip(strComp, losComp)]
    log_freq = [np.log10(f) for f in scaledFreq]
    log_strComp = [np.log10(s) for s in strComp]
    log_losComp = [np.log10(l) for l in losComp]

    if select == 0:
        y1 = strComp
        y2 = losComp
        y_label = "J', "+'J" /Pa'
        ylim = [np.max(y1)*(-0.1), np.max(y1)*1.1]
        label1 = 'Storage compliance (linear)'
        label2 = 'Loss compliance (linear)'
        c1 = 'r'
        c2 = 'b'
        a = 1
        legend_loc='upper right'
        savefile = './png/SLS1_complex_compliance_linear.png'

    if select == 1:
        y1 = log_strComp
        y2 = log_losComp
        y_label = "log[J', "+'J" /Pa]'
        label1 = 'Storage compliance (log)'
        label2 = 'Loss compliance (log)'
        c1 = 'r'
        c2 = 'b'
        a = 1
        legend_loc='upper right'
        savefile = './png/SLS1_complex_compliance_log.png'
        try:
            fitting = int(input('Selection (curve fit: 0, no curve fit: 1): '))
        except ValueError:
            fitting = 0
        if fitting == 1:
            pass
        if fitting == 0:
            fitFreqs = fitFreqs()
            fit_strComp, fit_losComp, fit_result1, fit_result2 = curveFit(freq, fitFreqs)

    if select == 2:
        y1 = losTan
        y2 = [0 for i in range(len(y1))]
        y_label = 'loss tangent /'
        ylim = [np.max(y1)*(-0.1), np.max(y1)*1.5]
        label1 = 'Loss tangent'
        label2 = ''
        c1 = 'g'
        c2 = 'b'
        a = 0
        legend_loc='upper right'
        savefile = './png/SLS1_loss_tangent.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title('SLS model I '+param_text)
    ax.set_xlabel(freq_label)
    ax.set_ylabel(y_label)
    ax.scatter(log_freq, y1, c=c1, label=label1)
    ax.scatter(log_freq, y2, c=c2, label=label2, alpha=a)

    if fitting == 1:
        ylim = [np.min(y1)-0.5, np.max(y1)+0.5]
    if fitting == 0:
        ax.plot(log_freq, fit_strComp, c='r', ls=':', label='fitted Storage compliance')
        ax.plot(log_freq, fit_losComp, c='b', ls=':', label='fitted Loss compliance')
        ylim = [np.min(y1)-0.5, np.max(y1)+0.5]
        fig.text(0.15, 0.30, fit_result1)
        fig.text(0.15, 0.25, fit_result2)

    ax.set_ylim(ylim[0], ylim[1])   
    ax.legend(loc=legend_loc)
    ax.grid()
    ax.set_axisbelow(True)
    fig.savefig(savefile)

    plt.show()