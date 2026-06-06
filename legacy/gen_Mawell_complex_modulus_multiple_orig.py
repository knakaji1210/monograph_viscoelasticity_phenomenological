# complex modulus of generalized Maxwell model (multiple freqs)
 
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib import gridspec as gs

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
    maxwellMod_list = []
    relaxTime_list = []
    try:
        infMod = float(input('Enter equilibrium modulus value (MPa) (default = 0.1 MPa): '))*10**6
    except ValueError:
        infMod = 10**5
    try:
        numComponent = int(input('Enter the number of Maxwell components (default = 1): '))
    except ValueError:
        numComponent = 1
    for i in range(numComponent):
        try:
            maxwellMod = float(input('Enter modulus value of Maxwell component (MPa) (default = 1 MPa): '))*10**6
        except ValueError:
            maxwellMod = 10**6
        maxwellMod_list.append(maxwellMod)
        try:
            viscosity = float(input('Enter viscosity value of Maxwell component (kPa s) (default = 100 kPa s): '))*10**3
        except ValueError:
            viscosity = 10**5
        relaxTime = viscosity/maxwellMod
        relaxTime_list.append(relaxTime)
    return numComponent, maxwellMod_list, infMod, relaxTime_list

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
    return freq

def fitFreqs():
    try:
        minFreq = float(input('Enter minimum frequency for fitting (storage) (default = 0.5): '))
    except ValueError:
        minFreq = 0.5
    try:
        maxFreq = float(input('Enter maximum frequency for fitting (storage) (default = 3): '))
    except ValueError:
         maxFreq = 3
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
    param,_ = curve_fit(loglogFit, log_freq[minFit:maxFit], log_strMod[minFit:maxFit])
    fit_strMod = fittedArray(log_freq, param)
    fit_result = "E' ∝ (freq*T)^({0:.2f})".format(param[0])
    return fit_strMod, fit_result

if __name__=='__main__':
    # calcul1ating complex modulus and loss tangent
    numComponent, maxwellMod_list, infMod, relaxTime_list = reqParams()
    param_text = '(E2 = {0:.1f} MPa, Num. Component = {1})'.format(infMod/10**6, numComponent)
    freq = freqAxis()
    freq_label = 'log[freq]'
    strMod = [infMod for s in range(len(freq))]
    losMod = [0 for l in range(len(freq))]
    for i in range(numComponent):
        strMod = [s1 + s2 for (s1, s2) in zip(strMod, complexMod(maxwellMod_list[i], 0, relaxTime_list[i], freq)[0])]
        losMod = [l1 + l2 for (l1, l2) in zip(losMod, complexMod(maxwellMod_list[i], 0, relaxTime_list[i], freq)[1])]
    losTan = [l/s for (s, l) in zip(strMod, losMod)]
    log_freq = [np.log10(f) for f in freq]
    log_strMod = [np.log10(s) for s in strMod]
    log_losMod = [np.log10(l) for l in losMod]
    strMod_c_list = []
    losMod_c_list = []
    losTan_c_list = []
    log_strMod_c_list = []
    log_losMod_c_list = []
    for i in range(numComponent):
        strMod_c = complexMod(maxwellMod_list[i], 0, relaxTime_list[i], freq)[0]
        losMod_c = complexMod(maxwellMod_list[i], 0, relaxTime_list[i], freq)[1]
        losTan_c = [l/s for (s, l) in zip(strMod_c, losMod_c)]
        log_strMod_c = [np.log10(s) for s in strMod_c]
        log_losMod_c = [np.log10(l) for l in losMod_c]
        strMod_c_list.append(strMod_c)
        losMod_c_list.append(losMod_c)
        losTan_c_list.append(losTan_c)
        log_strMod_c_list.append(log_strMod_c)
        log_losMod_c_list.append(log_losMod_c)
    fitting = -1

    try:
        select = int(input('Selection (complex modulus (linear): 0, complex modulus (log): 1, loss tangent: 2): '))
    except ValueError:
        select = 0

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
        savefile = './png/genMaxwell_complex_modulus_linear_m.png'

    if select == 1:
        y1 = log_strMod
        y2 = log_losMod
        y_label = "log[E', "+'E" /Pa]'
        label1 = 'Storage modulus (log)'
        label2 = 'Loss modulus (log)'
        c1 = 'r'
        c2 = 'b'
        a = 1
        savefile = './png/genMaxwell_complex_modulus_log_m.png'
        try:
            fitting = int(input('Selection (curve fit: 0, no curve fit: 1): '))
        except ValueError:
            fitting = 0
        if fitting == 1:
            pass
        if fitting == 0:
            fitFreqs = fitFreqs()
            fit_strMod, fit_result = curveFit(freq, fitFreqs)

    if select == 2:
        y1 = losTan
        y2 = [0 for i in range(len(y1))]
        y_label = 'loss tangent /'
#        ylim = [np.max(y1)*(-0.1), np.max(y1)*1.5]
        ylim = [-1, 50]
        label1 = 'Loss tangent'
        label2 = ''
        c1 = 'g'
        c2 = 'b'
        a = 0
        savefile = './png/genMaxwell_loss_tangent_m.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,6))
    spec = gs.GridSpec(ncols=4, nrows=2, width_ratios=[1,3,3,1],height_ratios=[1,5])
    ax = fig.add_subplot(spec[1, 0:4])
    ax.set_title('generalized Maxwell model '+param_text)
    ax.set_xlabel(freq_label)
    ax.set_ylabel(y_label)
    ax.scatter(log_freq, y1, c=c1, label=label1)
    ax.scatter(log_freq, y2, c=c2, label=label2, alpha=a)

    if select == 0:
        for i in range(numComponent):
            ax.plot(log_freq, strMod_c_list[i], c='r', linewidth=0.5)
            ax.plot(log_freq, losMod_c_list[i], c='b', linewidth=0.5)

    if select == 1:
        for i in range(numComponent):
            ax.plot(log_freq, log_strMod_c_list[i], c='r', linewidth=0.5)
            ax.plot(log_freq, log_losMod_c_list[i], c='b', linewidth=0.5)

    if select == 2:
        for i in range(numComponent):
            ax.plot(log_freq, losTan_c_list[i], c='g', linewidth=0.5)

    if fitting == 1:
        ylim = [np.min(y1)-0.5, np.max(y1)+0.5]
    if fitting == 0:
        ax.plot(log_freq, fit_strMod, c='b', ls=':', label='fitted Storage modulus')
        ylim = [np.min(y1)-0.5, np.max(y1)+0.5]
        fig.text(0.15, 0.45, fit_result)
        ax.set_xlabel(freq_label)

    ax.set_ylim(ylim[0], ylim[1])   
    ax.legend(loc='upper left') 
    ax.grid()
    ax.set_axisbelow(True)

    ax2 = fig.add_subplot(spec[0, 1:4])
    plt.axis('off')
    column_names=['{}'.format(i+1) for i in range(numComponent)]
    column_labels=[""," log(E1)", "log(1/T)"]
    log_maxwellMod_list = ['{0:.1f}'.format(np.log10(m)) for m in maxwellMod_list]
    log_invRelaxTime_list = ['{0:.1f}'.format(np.log10(1/t)) for t in relaxTime_list]
    tableData=[column_names, log_maxwellMod_list, log_invRelaxTime_list]
    ax2.table(cellText=tableData,rowLabels=column_labels,loc="center")

    fig.savefig(savefile)

    plt.show()