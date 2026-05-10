# complex compliance of Voigt model

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
    strComp = np.array(strComp)
    losComp = np.array(losComp)
    complexComp = [strComp, losComp]
    return complexComp

def reqParams():
    try:
        infMod = float(input('Enter equilibrium modulus value (MPa) (default = 1.0 MPa): '))*10**6
    except ValueError:
        infMod = 10**6
    try:
        viscosity = float(input('Enter viscosity value (kPa s) (default = 100.0 kPa s): '))*10**3
    except ValueError:
        viscosity = 10**5
    relaxTime = viscosity/infMod
    return infMod, relaxTime

def freqAxis(relaxTime):
    centerFreq = 1 / relaxTime
    freq = np.logspace(int(np.log10(centerFreq))-3, int(np.log10(centerFreq))+3, 49)
#    scaledFreq =([f*relaxTime for f in freq])    # 初期のバージョン
    scaledFreq = freq*relaxTime                  # freqはnp.ndarrayなのでこのように書ける
    freqAxis = [freq, scaledFreq]
    return freqAxis

def fitFreqs():
    try:
        minFreq1 = float(input('Enter minimum frequency for fitting (storage) (default = 1): '))
    except ValueError:
        minFreq1 = 1
    try:
        maxFreq1 = float(input('Enter maximum frequency for fitting (storage) (default = 3): '))
    except ValueError:
         maxFreq1 = 3
    try:
        minFreq2 = float(input('Enter minimum frequency for fitting (loss) (default = -3): '))
    except ValueError:
        minFreq2 = -3
    try:
        maxFreq2 = float(input('Enter maximum frequency for fitting (loss) (default = -2): '))
    except ValueError:
         maxFreq2 = -2
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
    paramStr,_ = curve_fit(loglogFit, log_freq[minFit1:maxFit1], log_strComp[minFit1:maxFit1])
    paramLos,_ = curve_fit(loglogFit, log_freq[minFit2:maxFit2], log_losComp[minFit2:maxFit2])
    fit_strComp = fittedArray(log_freq, paramStr)
    fit_losComp = fittedArray(log_freq, paramLos)
    fit_result1 = r"$J^\prime \propto (\omega\tau)^{{{0:.2f}}}$".format(paramStr[0])
    fit_result2 = r"$J^{{\prime\prime}} \propto (\omega\tau)^{{{0:.2f}}}$".format(paramLos[0])
    return fit_strComp, fit_losComp, fit_result1, fit_result2

if __name__=='__main__':
    # calcul1ating dynamic compliance and loss tangent
    infMod, relaxTime = reqParams()
    param_text = r'($E_{{\infty}}$ = {0:.1f} MPa, $\tau$ = {1:.1f} ms)'.format(infMod/10**6, relaxTime*10**3)
    freqAxis = freqAxis(relaxTime)
    fitting = -1
    try:
        select = int(input('Selection (complex compliance (linear): 0, complex compliance (log): 1, loss tangent: 2): '))
    except ValueError:
        select = 0

    freq = freqAxis[0]
    scaledFreq = freqAxis[1]
    freq_label = r'log($\omega\tau$)'
    strComp = complexComp(infMod, relaxTime, freq)[0]
    losComp = complexComp(infMod, relaxTime, freq)[1]
    losTan = np.array([l/s for (s, l) in zip(strComp, losComp)])
    log_freq = np.array([np.log10(f) for f in scaledFreq])
    log_strComp = np.array([np.log10(s) for s in strComp])
    log_losComp = np.array([np.log10(l) for l in losComp])

    if select == 0:
        y1 = strComp*10**6             # rescale to MPa
        y2 = losComp*10**6             # rescale to MPa
        y_label = r'$J^\prime$, $J^{{\prime\prime}}$ /MPa$^{{{-1}}}$'
        ylim = [np.max(y1)*(-0.1), np.max(y1)*1.1]
        label1 = r'$J^\prime$ (linear)'
        label2 = r'$J^{{\prime\prime}}$ (linear)'
        c1 = 'r'
        c2 = 'b'
        a = 1
        legend_loc='upper right'
        savefile = './png/Voigt_complex_compliance_linear.png'

    if select == 1:
        y1 = log_strComp
        y2 = log_losComp
        y_label = "log[J', "+'J" /Pa]'
        y_label = r'log($J^\prime$, $J^{{\prime\prime}}$ /Pa$^{{{-1}}}$)'
        label1 = r'$J^\prime$ (log)'
        label2 = r'$J^{{\prime\prime}}$ (log)'
        c1 = 'r'
        c2 = 'b'
        a = 1
        legend_loc='upper right'
        savefile = './png/Voigt_complex_compliance_log.png'
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
        y_label = 'tan $\delta$ /'
        ylim = [np.max(y1)*(-0.1), np.max(y1)*1.1]
        label1 = r'tan $\delta$'
        label2 = ''
        c1 = 'g'
        c2 = 'b'
        a = 0
        legend_loc='upper left'
        savefile = './png/Voigt_loss_tangent.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title('Voigt model '+param_text)
    ax.set_xlabel(freq_label)
    ax.set_ylabel(y_label)
    ax.scatter(log_freq, y1, c=c1, label=label1)
    ax.scatter(log_freq, y2, c=c2, label=label2, alpha=a)

    if fitting == 1:
        ylim = [np.min(y1)-1, np.max(y1)+1]
    if fitting == 0:
        ax.plot(log_freq, fit_strComp, c='r', ls=':', label=r'fitted $J^{\prime}$')
        ax.plot(log_freq, fit_losComp, c='b', ls=':', label=r'fitted $J^{{\prime\prime}}$')
        ylim = [np.min(y1)-1, np.max(y1)+1]
        fig.text(0.2, 0.30, fit_result1)
        fig.text(0.2, 0.25, fit_result2)

    ax.set_ylim(ylim[0], ylim[1])   
    ax.legend(loc=legend_loc)
    ax.grid()
    ax.set_axisbelow(True)
    fig.savefig(savefile, dpi=300)

    plt.show()