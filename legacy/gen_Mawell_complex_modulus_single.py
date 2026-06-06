# complex modulus of generalized Maxwell model (single freq)
 
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_complexMod(E1, E2, T, f):
    # E1: Maxwell spring component, E2: infMod
    numer = E1*T*f*(2j/2)
    denom = 1 + T*f*(2j/2)
    comMod = E2 + numer/denom
    strMod = comMod.real
    losMod = comMod.imag
    return strMod, losMod

def complexMod(E1, E2, T, freq):
    strMod = []
    losMod = []
    for f in freq:
        s, l = calc_complexMod(E1, E2, T, f)
        strMod.append(s)
        losMod.append(l)
    complexMod = [strMod, losMod]
    return complexMod

def reqParams():
    try:
        infMod = float(input('Enter equilibrium modulus value (MPa) (default = 0.1 MPa): '))*10**6
    except ValueError:
        infMod = 10**5
    try:
        maxwellMod = float(input('Enter modulus value of Maxwell component (MPa) (default = 1 MPa): '))*10**6
    except ValueError:
        maxwellMod = 10**6
    try:
        viscosity = float(input('Enter viscosity value of Maxwell component (kPa s) (default = 100 kPa s): '))*10**3
    except ValueError:
        viscosity = 10**5
    relaxTime = viscosity/maxwellMod
    return maxwellMod, infMod, relaxTime


def reqFreqs():
    try:
        minFreq = int(input('Enter minimum frequency in log scale (default = -3): '))
    except ValueError:
        minFreq = -3
    try:
        maxFreq = int(input('Enter maximum frequency in log scale (default = 6): '))
    except ValueError:
        maxFreq = 6
    intFreq = maxFreq - minFreq
    freqInfo = [minFreq, maxFreq, intFreq]
    return freqInfo

def freqAxis():
    freqInfo = reqFreqs()
    freq = np.logspace(freqInfo[0], freqInfo[1], freqInfo[2]*5+1)
    scaledFreq = [f*relaxTime for f in freq]
    freqAxis = [freq, scaledFreq]
    return freqAxis

def fitFreqs():
    try:
        minFreq = float(input('Enter minimum frequency for fitting (storage) (default = -3): '))
    except ValueError:
        minFreq = -3
    try:
        maxFreq = float(input('Enter maximum frequency for fitting (storage) (default = 6): '))
    except ValueError:
         maxFreq = 6
    fitFreqs = [minFreq, maxFreq]
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
    # for generalized model, scaledFreq was replaced to freq
    minFit = fitRegion(freq, fitFreqs[0], fitFreqs[1])[0]
    maxFit = fitRegion(freq, fitFreqs[0], fitFreqs[1])[1]
    param,_ = curve_fit(loglogFit, log_freq[minFit:maxFit], log_losMod[minFit:maxFit])
    fit_losMod = fittedArray(log_freq, param)
    fit_result = 'E" ∝ (freq*T)^({0:.2f})'.format(param[0])
    return fit_losMod, fit_result

if __name__=='__main__':
    # calcul1ating complex modulus and loss tangent
    maxwellMod, infMod, relaxTime = reqParams()
    param_text = '(E1 = {0:.1f} MPa, E2 = {1:.1f} MPa, T = {2:.1f} ms)'.format(maxwellMod/10**6, infMod/10**6, relaxTime*10**3)
    freqAxis = freqAxis()
    fitting = -1
    try:
        select = int(input('Selection (complex modulus (linear): 0, complex modulus (log): 1, loss tangent: 2): '))
    except ValueError:
        select = 0

    freq = freqAxis[0]
    scaledFreq = freqAxis[1]
#    freq_label = 'log[freq*T]'
    freq_label = 'log[freq]'
    strMod = complexMod(maxwellMod, infMod, relaxTime, freq)[0]
    losMod = complexMod(maxwellMod, infMod, relaxTime, freq)[1]
    losTan = [l/s for (s, l) in zip(strMod, losMod)]
#    log_freq = [np.log10(f) for f in scaledFreq]
    log_freq = [np.log10(f) for f in freq]
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
        savefile = './png/genMaxwell_complex_modulus_linear_s.png'

    if select == 1:
        y1 = log_strMod
        y2 = log_losMod
        y_label = "log[E', "+'E" /Pa]'
        label1 = 'Storage modulus (log)'
        label2 = 'Loss modulus (log)'
        c1 = 'r'
        c2 = 'b'
        a = 1
        savefile = './png/genMaxwell_complex_modulus_log_s.png'
        try:
            fitting = int(input('Selection (curve fit: 0, no curve fit: 1): '))
        except ValueError:
            fitting = 0
        if fitting == 1:
            pass
        if fitting == 0:
            fitFreqs = fitFreqs()
            fit_losMod, fit_result = curveFit(freq, fitFreqs)

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
        savefile = './png/genMaxwell_loss_tangent_s.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(111)
    ax.set_title('SLS model II '+param_text)
    ax.set_xlabel(freq_label)
    ax.set_ylabel(y_label)
    ax.scatter(log_freq, y1, c=c1, label=label1)
    ax.scatter(log_freq, y2, c=c2, label=label2, alpha=a)

    if fitting == 1:
        ylim = [np.min(y1)-0.5, np.max(y1)+0.5]
    if fitting == 0:
        ax.plot(log_freq, fit_losMod, c='b', ls=':', label='fitted Loss modulus')
        ylim = [np.min(y1)-0.5, np.max(y1)+0.5]
        fig.text(0.15, 0.45, fit_result)

    ax.set_ylim(ylim[0], ylim[1])   
    ax.legend(loc='upper left') 
    ax.grid()
    ax.set_axisbelow(True)
#    fig.savefig(savefile)

    plt.show()