# relaxation modulus of generalized Maxwell model (multiple freqs)

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# マクスウェル要素の緩和弾性率の計算式
def calcRelaxMod(E1, E2, tau, t):
    # E1: Maxwell spring component, E2: infMod
    relaxMod = E2 + E1*np.exp(-t/tau)
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
    E_list = []
    eta_list = []
    tau_list = []
    try:
        infMod = float(input('Enter equilibrium modulus value (MPa) (default = 0.1 MPa): '))*10**6
    except ValueError:
        infMod = 10**5
    try:
        numComp = int(input('Enter the number of Maxwell components (default = 1): '))
    except ValueError:
        numComp = 1
    for j in range(numComp):
        try:
            Ej = float(input('Enter modulus value of Maxwell component (MPa) (default = 1 MPa): '))*10**6
        except ValueError:
            Ej = 10**6
        E_list.append(Ej)
        try:
            etaj = float(input('Enter viscosity value of Maxwell component (kPa s) (default = 100 kPa s): '))*10**3
        except ValueError:
            etaj = 10**5
        eta_list.append(etaj)
        tauj = etaj/Ej
        tau_list.append(tauj)
    return numComp, infMod, E_list, eta_list, tau_list

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
    fit_result = "$E$($t$) ∝ $t^{{{0:.2f}}}$".format(param[0])
    return fit_relaxMod, fit_result

def calcRelaxSpectrum(x, y):
    dydx = np.gradient(y, x)
    relaxSpectrum = -10**y * dydx
    with np.errstate(divide='ignore'):
        log_relaxSpectrum = np.log10(relaxSpectrum)
    return relaxSpectrum, log_relaxSpectrum

if __name__=='__main__':
    # calcul1ating relaxation modulus
    select, time = timeAxisChoice()
    numComp, infMod, E_list, eta_list, tau_list = reqParams()
    insMod = infMod + np.sum(E_list)   # 瞬間弾性率=単独バネの（平衡）弾性率＋各マクスウェル要素の弾性率の和　式(6.12)
    param_text = r'($E_i$ = {0:.1f} MPa, $E_\infty$ = {1:.1f} MPa, {2} components)'.format(insMod/10**6, infMod/10**6, numComp)
    fitting = -1
    spectrum = -1

    # 緩和弾性率の計算
    relaxMod = infMod * np.ones(len(time))      # ゲタとして平衡弾性率のnp配列を用意
    for j in range(numComp):   # 各マクスウェル要素分を追加（ただし、平衡男性率は追加済みなのでcalcRelaxModの該当引数は0にする）
        relaxMod = relaxMod + calcRelaxMod(E_list[j], 0, tau_list[j], time)
    log_time = np.log10(time)
    with np.errstate(divide='ignore'):
        log_relaxMod = np.log10(relaxMod)
    # 各マクスウェル要素の緩和弾性率の計算（ただし全ての要素でゲタとして平衡弾性率を加える）
    relaxMod_comp = np.zeros((numComp,len(time)))
    for j in range(numComp):
        relaxMod_comp[j] = calcRelaxMod(E_list[j], infMod, tau_list[j], time)
    with np.errstate(divide='ignore'):
        log_relaxMod_comp = np.log10(relaxMod_comp)

    if select == 0:
        x = time
        y = log_relaxMod
        x_label = r'$t$ /s'
        y_label = r'log[$E$($t$) /Pa]'
        xlim = [0, np.max(x)]
        ylim = [np.min(y)-0.5, np.max(y)+0.5]
        label = r'Relaxation modulus (linear)'
        legend_loc = 'upper right'
        c = 'r'
        savefile = './png/gen_Maxwell_relaxMod_linear.png'

    if select == 1:
        x = log_time
        y = log_relaxMod
        ymin_l = np.max([np.min(y)-0.5,-0.5])
        y_orig = relaxMod
        x_label = r'log[$t$ /s]'
        y_label = r'log[$E$($t$) /Pa]'
        xlim = [np.min(x)-0.5, np.max(x)+0.5]
        ylim = [ymin_l, np.max(y)+0.5]
        label = r'Relaxation modulus (log)'
        legend_loc = 'upper right'
        c = 'r'
        savefile = './png/gen_Maxwell_relaxMod_log.png'
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
            relaxSpectrum, log_relaxSpectrum = calcRelaxSpectrum(x,y)

    # drawing graphs
    fig, ax = plt.subplots(figsize=(8, 6), tight_layout=True)
    ax.set_title('generalized Maxwell model '+param_text)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.scatter(x, y, c=c, label=label)

    for j in range(numComp):
        ax.plot(x, log_relaxMod_comp[j], c='b', linewidth=0.5)

    if fitting == 1:
        pass
    if fitting == 0:
        ax.plot(x, fit_relaxMod, c='b', ls=':', label='fitted Relaxation modulus')
        fig.text(0.5, 0.75, fit_result)
        ax.set_xlabel(x_label)

    if spectrum == 0:
        y = [rs/insMod for rs in relaxSpectrum]
        legend_loc='upper left'
        ax1 = ax.twinx()
        ax1.scatter(x, y, c='C0', ls=':', label='Relaxation Spectrum')
        ax1.set_ylim(-0.1*np.max(y),1.5*np.max(y))
        ax1.set_ylabel('Relaxation Spectrum / Eins')
        ax1.legend(loc='upper right')
        savefile = './png/genMaxwell_relaxation_modulus_log_spectrum.png'

    ax.set_xlim(xlim[0], xlim[1]) 
    ax.set_ylim(ylim[0], ylim[1])
    ax.legend(loc=legend_loc)
    ax.grid()
    ax.set_axisbelow(True)

    param_table = ax.table(
        cellText=[['{0:.1f}'.format(np.log10(E_list[j])) for j in range(numComp)],
                  ['{0:.1f}'.format(eta_list[j]/10**3) for j in range(numComp)],
                  ['{0:.2f}'.format(np.log10(tau_list[j])) for j in range(numComp)]],
        rowLabels=[r"log($E_j$ /Pa)", r"$\eta_j$ /kPa s", r"log($\tau_j$ /s)"],
        colLabels=['#{}'.format(j+1) for j in range(numComp)],
        loc='bottom',
        bbox=[0.0, -0.45, 1.0, 0.3]  # [x, y, 幅, 高さ] で微調整
    )
    param_table.set_fontsize(11)
    param_table.scale(1, 1.5) # セルの大きさを調整
    plt.subplots_adjust(bottom=0.2)

    fig.savefig(savefile, dpi=300)

    plt.show()