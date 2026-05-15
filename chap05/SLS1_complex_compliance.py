# complex compliance of SLS I model

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_complexComp(E, k, tau, af):
    numer = (k + tau*af*(2j/2))/E
    # E must be insMod
    denom = 1 + tau*af*(2j/2)
    comComp = numer/denom
    strComp = comComp.real
    losComp = -comComp.imag
    return strComp, losComp

def complexComp(E, k, tau, angFreq):
    strComp, losComp = calc_complexComp(E, k, tau, angFreq)
    complexComp = [strComp, losComp]
    return complexComp

def reqParams():
    try:
        E1 = float(input('modulus 1 [MPa] (default = 1.0 MPa): '))*10**6
    except ValueError:
        E1 = 10**6                  # [Pa] modulus
    try:
        E2 = float(input('modulus 2 [MPa] (default = 0.2 MPa): '))*10**6
    except ValueError:
        E2 = 2*10**5                # [Pa] modulus
    try:
        eta = float(input('viscosity [kPa s] (default = 500.0 kPa s): '))*10**3
    except ValueError:
        eta = 5*10**5               # [Pa s] viscosity

    insMod = E1                     # [Pa] instantaneous modulus
    infMod = E1*E2/(E1+E2)          # [Pa] equilibrium modulus
    k = insMod/infMod
    retardTime = eta/E2             # [s] retardation time
    return E1, E2, eta, insMod, infMod, k, retardTime


def freqAxis(retardTime):
    centerAngFreq = 1 / retardTime
    angFreq = np.logspace(int(np.log10(centerAngFreq))-2.5, int(np.log10(centerAngFreq))+1.5, 51)
    scaledAngFreq = angFreq*retardTime
    freqAxes = [angFreq, scaledAngFreq]
    return freqAxes

def fitAngFreqs():
    try:
        minAngFreq_s = float(input('Enter minimum frequency for fitting (storage) (default = 0.0): '))
    except ValueError:
        minAngFreq_s = 0.0
    try:
        maxAngFreq_s = float(input('Enter maximum frequency for fitting (storage) (default = 0.5): '))
    except ValueError:
         maxAngFreq_s = 0.5
    try:
        minAngFreq_l = float(input('Enter minimum frequency for fitting (loss) (default = 1.0): '))
    except ValueError:
        minAngFreq_l = 1.0
    try:
        maxAngFreq_l = float(input('Enter maximum frequency for fitting (loss) (default = 1.8): '))
    except ValueError:
         maxAngFreq_l = 1.8
    fitAngFreqs = [minAngFreq_s, maxAngFreq_s, minAngFreq_l, maxAngFreq_l]
    return fitAngFreqs

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

def curveFit(scaledAngFreq, fitAngFreqs):
    minFit1 = fitRegion(scaledAngFreq, fitAngFreqs[0], fitAngFreqs[1])[0]
    maxFit1 = fitRegion(scaledAngFreq, fitAngFreqs[0], fitAngFreqs[1])[1]
    minFit2 = fitRegion(scaledAngFreq, fitAngFreqs[2], fitAngFreqs[3])[0]
    maxFit2 = fitRegion(scaledAngFreq, fitAngFreqs[2], fitAngFreqs[3])[1]
    paramStr,_ = curve_fit(loglogFit, log_scaledAngFreq[minFit1:maxFit1], log_strComp[minFit1:maxFit1])
    paramLos,_ = curve_fit(loglogFit, log_scaledAngFreq[minFit2:maxFit2], log_losComp[minFit2:maxFit2])
    fit_strComp = fittedArray(log_scaledAngFreq, paramStr)
    fit_losComp = fittedArray(log_scaledAngFreq, paramLos)
    fit_result1 = r"$J^\prime \propto (\omega\tau)^{{{0:.2f}}}$".format(paramStr[0])
    fit_result2 = r"$J^{{\prime\prime}} \propto (\omega\tau)^{{{0:.2f}}}$".format(paramLos[0])
    return fit_strComp, fit_losComp, fit_result1, fit_result2

if __name__=='__main__':
    # calcul1ating complex compliance and loss tangent
    E1, E2, eta, insMod, infMod, k, retardTime = reqParams()
    param_text = r'($E_i$ = {0:.2f} MPa, $E_{{\infty}}$ = {1:.2f} MPa, $\tau$ = {2:.1f} s)'.format(insMod/10**6, infMod/10**6, retardTime)
    freqAxes = freqAxis(retardTime)
    fitting = -1
    try:
        select = int(input('Selection (complex compliance (linear): 0, complex compliance (log): 1, loss tangent: 2): '))
    except ValueError:
        select = 0

    angFreq = freqAxes[0]
    scaledAngFreq = freqAxes[1]
    scaledAngFreq_label = r'log($\omega\tau$)'
    strComp = complexComp(insMod, k, retardTime, angFreq)[0]
    losComp = complexComp(insMod, k, retardTime, angFreq)[1]
    losTan = losComp / strComp
    log_scaledAngFreq = np.log10(scaledAngFreq)
    log_strComp = np.log10(strComp)
    log_losComp = np.log10(losComp)

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
        savefile = './png/SLS1_complex_compliance_linear.png'

    if select == 1:
        y1 = log_strComp
        y2 = log_losComp
        y_label = r'log($J^\prime$, $J^{{\prime\prime}}$ /Pa$^{{{-1}}}$)'
        ylim = [np.min(y1)-3, np.max(y1)+1]
        label1 = r'$J^\prime$ (log)'
        label2 = r'$J^{{\prime\prime}}$ (log)'
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
            fitAngFreqs = fitAngFreqs()
            fit_strComp, fit_losComp, fit_result1, fit_result2 = curveFit(scaledAngFreq, fitAngFreqs)

    if select == 2:
        y1 = losTan
        y2 = [0 for i in range(len(y1))]
        y_label = 'tan $\delta$ /'
        ylim = [np.max(y1)*(-0.1), np.max(y1)*1.5]
        label1 = r'tan $\delta$'
        label2 = ''
        c1 = 'g'
        c2 = 'b'
        a = 0
        legend_loc='upper right'
        savefile = './png/SLS1_loss_tangent.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title('SLS I model '+param_text)
    ax.set_xlabel(scaledAngFreq_label)
    ax.set_ylabel(y_label)
    ax.scatter(log_scaledAngFreq, y1, c=c1, label=label1)
    ax.scatter(log_scaledAngFreq, y2, c=c2, label=label2, alpha=a)

    if fitting == 1:
        pass
    if fitting == 0:
        ax.plot(log_scaledAngFreq, fit_strComp, c='r', ls=':', label=r'fitted $J^{\prime}$')
        ax.plot(log_scaledAngFreq, fit_losComp, c='b', ls=':', label=r'fitted $J^{{\prime\prime}}$')
        fig.text(0.2, 0.30, fit_result1)
        fig.text(0.2, 0.25, fit_result2)

    ax.set_ylim(ylim[0], ylim[1])   
    ax.legend(loc=legend_loc)
    ax.grid()
    ax.set_axisbelow(True)
    fig.savefig(savefile)

    plt.show()