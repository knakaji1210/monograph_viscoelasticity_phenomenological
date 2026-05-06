# complex modulus of Maxwell model

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_complexMod(E, T, f):
    numer = E*T*f*(2j/2)
    denom = 1 + T*f*(2j/2)
    comMod = numer/denom
    strMod = comMod.real
    losMod = comMod.imag
    return strMod, losMod

def complexMod(E, T, freq):
    strMod = []
    losMod = []
    for f in freq:
        s, l = calc_complexMod(E, T, f)
        strMod.append(s)
        losMod.append(l)
    strMod = np.array(strMod)
    losMod = np.array(losMod)
    dynamicMod = [strMod, losMod]
    return dynamicMod

def reqParams():
    try:
        insMod = float(input('Enter instantaneous modulus value (MPa) (default = 1.0 MPa): '))*10**6
    except ValueError:
        insMod = 10**6
    try:
        viscosity = float(input('Enter viscosity value (kPa s) (default = 100.0 kPa s): '))*10**3
    except ValueError:
        viscosity = 10**5
    relaxTime = viscosity/insMod
    return insMod, relaxTime

def freqAxis(relaxTime):
    centerFreq = 1 / relaxTime
    freq = np.logspace(int(np.log10(centerFreq))-3, int(np.log10(centerFreq))+3, 49)
#    scaledFreq = [f*relaxTime for f in freq]    # 初期のバージョン
    scaledFreq = freq*relaxTime                  # freqはnp.ndarrayなのでこのように書ける
    freqAxis = [freq, scaledFreq]
    return freqAxis

def fitFreqs():
    try:
        minFreq1 = float(input('Enter minimum frequency for fitting (storage) (default = -3): '))
    except ValueError:
        minFreq1 = -3
    try:
        maxFreq1 = float(input('Enter maximum frequency for fitting (storage) (default = -1): '))
    except ValueError:
         maxFreq1 = -1
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
    fitted_array = np.array(fitted_array)
    return fitted_array

def curveFit(freq, fitFreqs):
    minFit1 = fitRegion(scaledFreq, fitFreqs[0], fitFreqs[1])[0]
    maxFit1 = fitRegion(scaledFreq, fitFreqs[0], fitFreqs[1])[1]
    minFit2 = fitRegion(scaledFreq, fitFreqs[2], fitFreqs[3])[0]
    maxFit2 = fitRegion(scaledFreq, fitFreqs[2], fitFreqs[3])[1]
    paramStr,_ = curve_fit(loglogFit, log_freq[minFit1:maxFit1], log_strMod[minFit1:maxFit1])
    paramLos,_ = curve_fit(loglogFit, log_freq[minFit2:maxFit2], log_losMod[minFit2:maxFit2])
    fit_strMod = fittedArray(log_freq, paramStr)
    fit_losMod = fittedArray(log_freq, paramLos)
    fit_result1 = r"$E^\prime \propto (\omega\tau)^{{{0:.2f}}}$".format(paramStr[0])
    fit_result2 = r"$E^{{\prime\prime}} \propto (\omega\tau)^{{{0:.2f}}}$".format(paramLos[0])
    return fit_strMod, fit_losMod, fit_result1, fit_result2

if __name__=='__main__':
    # calcul1ating dynamic Modulus and loss tangent
    insMod, relaxTime = reqParams()
    param_text = r'($E_i$ = {0:.1f} MPa, $\tau$ = {1:.1f} ms)'.format(insMod/10**6, relaxTime*10**3)
    freqAxis = freqAxis(relaxTime)
    fitting = -1
    try:
        select = int(input('Selection (complex modulus (linear): 0, complex modulus (log): 1, loss tangent: 2): '))
    except ValueError:
        select = 0

    freq = freqAxis[0]
    scaledFreq = freqAxis[1]
    freq_label = r'log($\omega\tau$)'
    strMod = complexMod(insMod, relaxTime, freq)[0]
    losMod = complexMod(insMod, relaxTime, freq)[1]
    losTan = np.array([l/s for (s, l) in zip(strMod, losMod)])
    log_freq = np.array([np.log10(f) for f in scaledFreq])
    log_strMod = np.array([np.log10(s) for s in strMod])
    log_losMod = np.array([np.log10(l) for l in losMod])

    if select == 0:
        y1 = strMod/10**6             # rescale to MPa
        y2 = losMod/10**6             # rescale to MPa
        y_label = r'$E^\prime$, $E^{{\prime\prime}}$ /MPa'
        ylim = [np.max(y1)*(-0.1), np.max(y1)*1.1]
        label1 = r'$E^\prime$ (linear)'
        label2 = r'$E^{{\prime\prime}}$ (linear)'
        c1 = 'r'
        c2 = 'b'
        a = 1
        legend_loc='upper left'
        savefile = './png/Maxwell_complex_modulus_linear.png'

    if select == 1:
        y1 = log_strMod
        y2 = log_losMod
        y_label = r'log($E^\prime$, $E^{{\prime\prime}}$ /Pa)'
        label1 = r'$E^\prime$ (log)'
        label2 = r'$E^{{\prime\prime}}$ (log)'
        c1 = 'r'
        c2 = 'b'
        a = 1
        legend_loc='upper left'
        savefile = './png/Maxwell_complex_modulus_log.png'
        try:
            fitting = int(input('Selection (curve fit: 0, no curve fit: 1): '))
        except ValueError:
            fitting = 0
        if fitting == 1:
            pass
        if fitting == 0:
            fitFreqs = fitFreqs()
            fit_strMod, fit_losMod, fit_result1, fit_result2 = curveFit(freq, fitFreqs)

    if select == 2:
        y1 = losTan
        y2 = np.zeros(len(y1))
        y_label = 'tan $\delta$ /'
        ylim = [np.max(y1)*(-0.1), np.max(y1)*1.1]
        label1 = r'tan $\delta$'
        label2 = ''
        c1 = 'g'
        c2 = 'b'
        a = 0
        legend_loc='upper right'
        savefile = './png/Maxwell_loss_tangent.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title('Maxwell model '+param_text)
    ax.set_xlabel(freq_label)
    ax.set_ylabel(y_label)
    ax.scatter(log_freq, y1, c=c1, label=label1)
    ax.scatter(log_freq, y2, c=c2, label=label2, alpha=a)

    if fitting == 1:
        ylim = [np.min(y1)-1, np.max(y1)+1]
    if fitting == 0:
        ax.plot(log_freq, fit_strMod, c='r', ls=':', label=r'fitted $E^{\prime}$')
        ax.plot(log_freq, fit_losMod, c='b', ls=':', label=r'fitted $E^{{\prime\prime}}$')
        ylim = [np.min(y1)-1, np.max(y1)+1]
        fig.text(0.7, 0.30, fit_result1)
        fig.text(0.7, 0.25, fit_result2)

    ax.set_ylim(ylim[0], ylim[1])   
    ax.legend(loc=legend_loc)
    ax.grid()
    ax.set_axisbelow(True)
    fig.savefig(savefile, dpi=300)

    plt.show()