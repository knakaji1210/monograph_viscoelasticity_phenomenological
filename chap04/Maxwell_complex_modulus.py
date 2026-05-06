# complex modulus of Maxwell model

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def reqParams():
    try:
        insMod = float(input('Enter instantaneous modulus value (MPa) (default = 0.4 MPa): '))*10**6
    except ValueError:
        insMod = 4*10**5
    try:
        viscosity = float(input('Enter viscosity value (kPa s) (default = 100.0 kPa s): '))*10**3
    except ValueError:
        viscosity = 10**5
    relaxTime = viscosity/insMod
    return insMod, relaxTime

def freqAxes(relaxTime):
    centerAngFreq = 1 / relaxTime
    angFreq = np.logspace(int(np.log10(centerAngFreq))-1.5, int(np.log10(centerAngFreq))+2.5, 51)
    scaledAngFreq = angFreq*relaxTime
    freqAxes = [angFreq, scaledAngFreq]
    return freqAxes

def calc_complexMod(E, tau, af):
    numer = E*af*tau*(2j/2)
    denom = 1 + af*tau*(2j/2)
    comMod = numer/denom
    strMod = comMod.real
    losMod = comMod.imag
    return strMod, losMod

def complexMod(E, tau, angFreq):
    strMod, losMod = calc_complexMod(E, tau, angFreq)
    dynamicMod = [strMod, losMod]
    return dynamicMod

def fitAngFreqs():
    try:
        minAngFreq_s = float(input('Enter minimum frequency for fitting (storage) (default = -2.0): '))
    except ValueError:
        minAngFreq_s = -2.0
    try:
        maxAngFreq_s = float(input('Enter maximum frequency for fitting (storage) (default = -0.8): '))
    except ValueError:
        maxAngFreq_s = -0.8
    try:
        minAngFreq_l = float(input('Enter minimum frequency for fitting (loss) (default = -2.0): '))
    except ValueError:
        minAngFreq_l = -2.0
    try:
        maxAngFreq_l = float(input('Enter maximum frequency for fitting (loss) (default = -0.8): '))
    except ValueError:
         maxAngFreq_l = -0.8
    fitAngFreqs = [minAngFreq_s, maxAngFreq_s, minAngFreq_l, maxAngFreq_l]
    return fitAngFreqs

def getNearestIdx(list, num):
    idx = np.abs(np.asarray(list) - num).argmin()
    return idx

def fitRegion(angFreq, minNum, maxNum):
    minFit = getNearestIdx(angFreq, 10**minNum)
    maxFit = getNearestIdx(angFreq, 10**maxNum)
    fitRegion = [minFit, maxFit]
    return fitRegion

def loglogFit(x, a, b):
    return  a*x + b

def fittedArray(x_array, param):
    fitted_array = [loglogFit(num, param[0], param[1]) for num in x_array]
    fitted_array = np.array(fitted_array)
    return fitted_array

def curveFit(scaledAngFreq, fitAngFreqs):
    minFit1 = fitRegion(scaledAngFreq, fitAngFreqs[0], fitAngFreqs[1])[0]
    maxFit1 = fitRegion(scaledAngFreq, fitAngFreqs[0], fitAngFreqs[1])[1]
    minFit2 = fitRegion(scaledAngFreq, fitAngFreqs[2], fitAngFreqs[3])[0]
    maxFit2 = fitRegion(scaledAngFreq, fitAngFreqs[2], fitAngFreqs[3])[1]
    paramStr,_ = curve_fit(loglogFit, log_scaledAngFreq[minFit1:maxFit1], log_strMod[minFit1:maxFit1])
    paramLos,_ = curve_fit(loglogFit, log_scaledAngFreq[minFit2:maxFit2], log_losMod[minFit2:maxFit2])
    fit_strMod = fittedArray(log_scaledAngFreq, paramStr)
    fit_losMod = fittedArray(log_scaledAngFreq, paramLos)
    fit_result1 = r"$E^\prime \propto (\omega\tau)^{{{0:.2f}}}$".format(paramStr[0])
    fit_result2 = r"$E^{{\prime\prime}} \propto (\omega\tau)^{{{0:.2f}}}$".format(paramLos[0])
    return fit_strMod, fit_losMod, fit_result1, fit_result2

if __name__=='__main__':
    # calcul1ating dynamic Modulus and loss tangent
    insMod, relaxTime = reqParams()
    param_text = r'($E_i$ = {0:.1f} MPa, $\tau$ = {1:.2f} s)'.format(insMod/10**6, relaxTime)
    freqAxes = freqAxes(relaxTime)
    fitting = -1    # initial value

    try:
        select = int(input('Selection (complex modulus (linear): 0, complex modulus (log): 1, loss tangent: 2): '))
    except ValueError:
        select = 0

    angFreq = freqAxes[0]
    scaledAngFreq = freqAxes[1]
    scaledAngFreq_label = r'log($\omega\tau$)'
    strMod = complexMod(insMod, relaxTime, angFreq)[0]
    losMod = complexMod(insMod, relaxTime, angFreq)[1]
    losTan = losMod / strMod
    log_scaledAngFreq = np.log10(scaledAngFreq)
    log_strMod = np.log10(strMod)
    log_losMod = np.log10(losMod)

    if select == 0:
        y1 = strMod/10**6             # rescale to MPa
        y2 = losMod/10**6             # rescale to MPa
        y_label = r'$E^\prime$, $E^{{\prime\prime}}$ /MPa'
        ylim = [np.max(y1)*(-0.2), np.max(y1)*1.2]
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
        ylim = [np.min(y1)-2, np.max(y1)+2]
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
            fitAngFreqs = fitAngFreqs()
            fit_strMod, fit_losMod, fit_result1, fit_result2 = curveFit(scaledAngFreq, fitAngFreqs)

    if select == 2:
        y1 = losTan
        y2 = np.zeros(len(y1))
        y_label = 'tan $\delta$ /'
        ylim = [np.max(y1)*(-0.2), np.max(y1)*1.2]
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
    ax.set_xlabel(scaledAngFreq_label)
    ax.set_ylabel(y_label)
    ax.scatter(log_scaledAngFreq, y1, c=c1, label=label1)
    ax.scatter(log_scaledAngFreq, y2, c=c2, label=label2, alpha=a)

    if fitting == 1:
        pass
    if fitting == 0:
        ax.plot(log_scaledAngFreq, fit_strMod, c='r', lw=1, ls=':', label=r'fitted $E^{\prime}$')
        ax.plot(log_scaledAngFreq, fit_losMod, c='b', lw=1, ls=':', label=r'fitted $E^{{\prime\prime}}$')
        ylim = [np.min(y1)-2, np.max(y1)+2]
        fig.text(0.7, 0.30, fit_result1)
        fig.text(0.7, 0.25, fit_result2)

    ax.set_ylim(ylim[0], ylim[1])   
    ax.legend(loc=legend_loc)
    ax.grid()
    ax.set_axisbelow(True)
    fig.savefig(savefile, dpi=300)

    plt.show()