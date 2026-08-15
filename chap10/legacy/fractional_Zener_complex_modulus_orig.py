# complex modulus of fractional Maxwell model
 
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_complexMod(E, k, T, f, nu):
    numer = E*(1/k + (T*f*(2j/2))**nu)
    # E must be insMod
    denom = 1 + (T*f*(2j/2))**nu
    comMod = numer/denom
    strMod = comMod.real
    losMod = comMod.imag
    return strMod, losMod

def complexMod(E, k, T, freq, nu):
    strMod = []
    losMod = []
    for f in freq:
        s, l = calc_complexMod(E, k, T, f, nu)
        strMod.append(s)
        losMod.append(l)
    dynamicMod = [strMod, losMod]
    return dynamicMod

def reqParams():
    try:
        insMod = float(input('Enter instantaneous modulus value (GPa) (default = 1 GPa): '))*10**9
    except ValueError:
        insMod = 10**9
    try:
        infMod = float(input('Enter equilibrium modulus value (MPa) (default = 1 MPa): '))*10**6
    except ValueError:
        infMod = 10**6
    try:
        modulus = float(input('Enter modulus value of spring-pot (GPa) (default = 1 GPa): '))*10**9
    except ValueError:
        modulus = 10**9
    try:
        viscosity = float(input('Enter viscosity value of spring-pot (kPa s) (default = 1000 kPa s): '))*10**3
    except ValueError:
        viscosity = 10**6
    try:
        nu = float(input('Enter fractional power (0 < nu < 1) (default = 0.5): '))
    except ValueError:
        nu = 0.5
    k = insMod/infMod
    kappa = modulus / (insMod - infMod)
    relaxTime = kappa**(1/nu)*viscosity/modulus
    return insMod, infMod, k, modulus, relaxTime, nu

def freqAxis(relaxTime):
    centerFreq = 1 / relaxTime
    freq = np.logspace(int(np.log10(centerFreq))-9, int(np.log10(centerFreq))+5, 57)
    scaledFreq = [f*relaxTime for f in freq]
    freqAxis = [freq, scaledFreq]
    return freqAxis

def fitFreqs():
    try:
        minFreq1 = float(input('Enter minimum frequency for fitting (storage) (default = -4): '))
    except ValueError:
        minFreq1 = -4
    try:
        maxFreq1 = float(input('Enter maximum frequency for fitting (storage) (default = -1): '))
    except ValueError:
         maxFreq1 = -1
    try:
        minFreq2 = float(input('Enter minimum frequency for fitting (loss) (default = -8): '))
    except ValueError:
        minFreq2 = -8
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
    fit_result1 = "E' ∝ (freq*tau)^({0:.2f})".format(paramStr[0])
    fit_result2 = 'E" ∝ (freq*tau)^({0:.2f})'.format(paramLos[0])
    return fit_strMod, fit_losMod, fit_result1, fit_result2

if __name__=='__main__':
    # calcul1ating complex modulus and loss tangent
    insMod, infMod, k, modulus, relaxTime, nu = reqParams()
    print(relaxTime)
    param_text = '(Ei = {0:.1f} MPa, Einf = {1:.1f} MPa, E = {2:.1f} MPa, \n tau = {3:.1f} ms, nu = {4:.2f})'.format(insMod/10**6, infMod/10**6, modulus/10**6, relaxTime*10**3, nu)
    freqAxis = freqAxis(relaxTime)
    fitting = -1
    try:
        select = int(input('Selection (complex modulus (linear): 0, complex modulus (log): 1, loss tangent: 2): '))
    except ValueError:
        select = 0

    freq = freqAxis[0]
    scaledFreq = freqAxis[1]
    freq_label = 'log[freq*tau]'
    strMod = complexMod(insMod, k, relaxTime, freq, nu)[0]
    losMod = complexMod(insMod, k, relaxTime, freq, nu)[1]
    losTan = [l/s for (s, l) in zip(strMod, losMod)]
    log_freq = [np.log10(f) for f in scaledFreq]
    log_strMod = [np.log10(s) for s in strMod]
    log_losMod = [np.log10(l) for l in losMod]

    if select == 0:
        y1 = strMod
        y2 = losMod
        y_label = "E', "+'E" /Pa'
        ylim = [np.max(y1)*(-0.1), np.max(y1)*1.1]
        label1 = 'Storage modulus (linear)'
        label2 = 'Loss modulus (linear)'
        c1 = 'r'
        c2 = 'b'
        a = 1
        legend_loc='upper left'
        savefile = './png/fractional_Zener_complex_modulus_linear.png'

    if select == 1:
        y1 = log_strMod
        y2 = log_losMod
        y_label = "log[E', "+'E" /Pa]'
        label1 = 'Storage modulus (log)'
        label2 = 'Loss modulus (log)'
        c1 = 'r'
        c2 = 'b'
        a = 1
        legend_loc='upper left'
        savefile = './png/fractional_Zener_complex_modulus_log.png'
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
        y2 = [0 for i in range(len(y1))]
        y_label = 'loss tangent /'
        ylim = [np.max(y1)*(-0.1), np.max(y1)*1.5]
        label1 = 'Loss tangent'
        label2 = ''
        c1 = 'g'
        c2 = 'b'
        a = 0
        legend_loc='upper left'
        savefile = './png/fractional_Zener_loss_tangent.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title('fractional Zener model '+param_text)
    ax.set_xlabel(freq_label)
    ax.set_ylabel(y_label)
    ax.scatter(log_freq, y1, c=c1, label=label1)
    ax.scatter(log_freq, y2, c=c2, label=label2, alpha=a)

    if fitting == 1:
        ylim = [np.min(np.minimum(y1, y2))-0.5, np.max(np.maximum(y1, y2))+0.5]
    if fitting == 0:
        ax.plot(log_freq, fit_strMod, c='r', ls=':', label='fitted Storage modulus')
        ax.plot(log_freq, fit_losMod, c='b', ls=':', label='fitted Loss modulus')
        ylim = [np.min(np.minimum(y1, y2))-0.5, np.max(np.maximum(y1, y2))+0.5]
        fig.text(0.7, 0.35, fit_result1)
        fig.text(0.7, 0.30, fit_result2)

    ax.set_ylim(ylim[0], ylim[1])   
    ax.legend(loc=legend_loc)
    ax.grid()
    ax.set_axisbelow(True)
    fig.savefig(savefile)

    plt.show()