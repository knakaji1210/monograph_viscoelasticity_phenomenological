# relaxation modulus of generalized Maxwell model (multiple freqs)

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib import gridspec as gs

def calc_relaxMod(E1, E2, T, t):
    # E1: Maxwell spring component, E2: infMod
    relaxMod = E2 + E1*np.exp(-t/T)
    return relaxMod

def relaxModFunc(E1, E2, T, time):
    relaxMod = []
    for t in time:
        c = calc_relaxMod(E1, E2, T, t)
        relaxMod.append(c)
    return relaxMod

def reqTimes():
    try:
        minTime = int(input('Enter minimum time in log scale (default = -7): '))
    except ValueError:
        minTime = -7
    try:
        maxTime = int(input('Enter maximum time in log scale (default = 1): '))
    except ValueError:
        maxTime = 1
    intTime = maxTime - minTime
    timeInfo = [minTime, maxTime, intTime]
    return timeInfo

def timeAxis():
    timeInfo = reqTimes()
    linearTime = np.linspace(10**timeInfo[0], 10**timeInfo[1], timeInfo[2]*30+1)
    logTime = np.logspace(timeInfo[0], timeInfo[1], timeInfo[2]*10+1)
    timeAxis = [linearTime, logTime]
    return timeAxis

def timeAxisChoice():
    try:
        select = int(input('Selection (relaxation modulus (linear): 0, relaxation modulus (log): 1): '))
    except ValueError:
        select = 0
    if select == 0:
        time = timeAxis()[0]
    if select == 1:
        time = timeAxis()[1]
    return select, time

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

def fitTimes():
    try:
        minTime = float(input('Enter minimum time for fitting (default = 0.5): '))
    except ValueError:
        minTime = 0.5
    try:
        maxTime = float(input('Enter maximum fime for fitting (default = 3): '))
    except ValueError:
         maxTime = 3
    fitTimes = [minTime, maxTime]
    return fitTimes

def getNearestIdx(list, num):
    idx = np.abs(np.asarray(list) - num).argmin()
    return idx

def fitRegion(time, minNum, maxNum):
    minFit = getNearestIdx(time, 10**minNum)
    maxFit = getNearestIdx(time, 10**maxNum)
    fitRegion = [minFit, maxFit]
    return fitRegion

def loglogFit(x, a, b):
    return  a*x + b

def fittedArray(x_array, param):
    fitted_array = [loglogFit(num, param[0], param[1]) for num in x_array]
    return fitted_array

def curveFit(time, fitFreqs):
    minFit = fitRegion(time, fitTimes[0], fitTimes[1])[0]
    maxFit = fitRegion(time, fitTimes[0], fitTimes[1])[1]
    param,_ = curve_fit(loglogFit, log_time[minFit:maxFit], log_relaxMod[minFit:maxFit])
    fit_relaxMod = fittedArray(log_time, param)
    fit_result = "E(t) ∝ t^({0:.2f})".format(param[0])
    return fit_relaxMod, fit_result

def relaxSpectrumFunc(x, y, y_orig):
    y_nd = []
    for i in range(len(x)-1):
        if i == 0:
            pass
        else:
            nd = (y[i+1]-y[i-1])/(x[i+1]-x[i-1])
            y_nd.append(nd)   
    x_cropped = np.delete(x,-1)
    x_cropped = np.delete(x_cropped,0)
    y_cropped = np.delete(y,-1)
    y_cropped = np.delete(y_cropped,0)
    y_orig_cropped = np.delete(y_orig,-1)
    y_orig_cropped = np.delete(y_orig_cropped,0)
    relaxSpectrum = [-y*nd for (y, nd) in zip(y_orig_cropped, y_nd)]
    rs = [-y*nd for (y, nd) in zip(y_orig_cropped, y_nd)]
#    relaxSpectrum = np.nan_to_num(rs, nan=1e-100)
    with np.errstate(divide='ignore'):
        log_relaxSpectrum = [np.log10(r) for r in relaxSpectrum]
    return x_cropped, relaxSpectrum, log_relaxSpectrum

if __name__=='__main__':
    # calcul1ating relaxation modulus
    select, time = timeAxisChoice()
    numComponent, maxwellMod_list, infMod, relaxTime_list = reqParams()
    insMod = infMod + np.sum(maxwellMod_list)
    param_text = '(Einf = {0:.1f} MPa, Eins = {1:.1f} MPa, Num. Component = {2})'.format(infMod/10**6, insMod/10**6, numComponent)
    fitting = -1
    spectrum = -1

    relaxMod = [infMod for r in range(len(time))]
    for i in range(numComponent):
        relaxMod = [r1 + r2 for (r1, r2) in zip(relaxMod, relaxModFunc(maxwellMod_list[i], 0, relaxTime_list[i], time))]
    log_time = [np.log10(t) for t in time]
    with np.errstate(divide='ignore'):
        log_relaxMod = [np.log10(r) for r in relaxMod]
    relaxMod_c_list = []
    log_relaxMod_c_list = []
    for i in range(numComponent):
        relaxMod_c = relaxModFunc(maxwellMod_list[i], infMod, relaxTime_list[i], time)
        with np.errstate(divide='ignore'):
            log_relaxMod_c = [np.log10(r) for r in relaxMod_c]
        relaxMod_c_list.append(relaxMod_c)
        log_relaxMod_c_list.append(log_relaxMod_c)

    if select == 0:
        x = time
        y = log_relaxMod
        x_label = 't /s'
        y_label = 'log[E(t) /Pa]'
        xlim = [0, np.max(x)]
        ylim = [np.min(y)-0.5, np.max(y)+0.5]
        label = 'Relaxation modulus (linear)'
        legend_loc = 'upper right'
        c = 'r'
        savefile = './png/genMaxwell_relaxation_modulus_linear_m.png'

    if select == 1:
        x = log_time
        y = log_relaxMod
        ymin_l = np.max([np.min(y)-0.5,-0.5])
        y_orig = relaxMod
        x_label = 'log[t /s]'
        y_label = 'log[E(t) /Pa]'
        xlim = [np.min(x)-0.5, np.max(x)+0.5]
        ylim = [ymin_l, np.max(y)+0.5]
        label = 'Relaxation modulus (log)'
        legend_loc = 'upper right'
        c = 'r'
        savefile = './png/genMaxwell_relaxation_modulus_log_m.png'
        try:
            fitting = int(input('Selection (curve fit: 0, no curve fit: 1): '))
        except ValueError:
            fitting = 1
        if fitting == 1:
            pass
        if fitting == 0:
            fitTimes = fitTimes()
            fit_relaxMod, fit_result = curveFit(time, fitTimes)
        try:
            spectrum = int(input('Selection (with spectrum: 0, without spectrum: 1): '))
        except ValueError:
            spectrum = 1
        if spectrum == 1:
            pass
        if spectrum == 0:
            x_cropped, relaxSpectrum, log_relaxSpectrum = relaxSpectrumFunc(x,y,y_orig)
        savefile = './png/genMaxwell_relaxation_modulus_log_m.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,6), tight_layout=True)
    spec = gs.GridSpec(ncols=4, nrows=2, width_ratios=[1,3,3,1],height_ratios=[1,5])
    ax = fig.add_subplot(spec[1, 0:4])
    ax.set_title('generalized Maxwell model '+param_text)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.scatter(x, y, c=c, label=label)

    for i in range(numComponent):
        ax.plot(x, log_relaxMod_c_list[i], c='b', linewidth=0.5)

    if fitting == 1:
        pass
    if fitting == 0:
        ax.plot(x, fit_relaxMod, c='b', ls=':', label='fitted Relaxation modulus')
        fig.text(0.15, 0.20, fit_result)
        ax.set_xlabel(x_label)

    if spectrum == 0:
        y = [rs/insMod for rs in relaxSpectrum]
        legend_loc='upper left'
        ax1 = ax.twinx()
        ax1.scatter(x_cropped, y, c='C0', ls=':', label='Relaxation Spectrum')
        ax1.set_ylim(-0.1*np.max(y),1.5*np.max(y))
        ax1.set_ylabel('Relaxation Spectrum / Eins')
        ax1.legend(loc='upper right')
        savefile = './png/genMaxwell_relaxation_modulus_log_spectrum.png'

    ax.set_xlim(xlim[0], xlim[1]) 
    ax.set_ylim(ylim[0], ylim[1])
    ax.legend(loc=legend_loc)
    ax.grid()
    ax.set_axisbelow(True)

    ax2 = fig.add_subplot(spec[0, 1:4])
    plt.axis('off')
    column_names=['{}'.format(i+1) for i in range(numComponent)]
    column_labels=[""," log(E1)", "log(tau)"]
    log_maxwellMod_list = ['{0:.1f}'.format(np.log10(m)) for m in maxwellMod_list]
    log_relaxTime_list = ['{0:.1f}'.format(np.log10(t)) for t in relaxTime_list]
    tableData=[column_names, log_maxwellMod_list, log_relaxTime_list]
    ax2.table(cellText=tableData,rowLabels=column_labels,loc="center")

    fig.savefig(savefile)

    plt.show()