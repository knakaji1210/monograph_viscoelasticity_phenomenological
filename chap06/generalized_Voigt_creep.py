# creep compliance of generalized Voigt model (multiple freqs)

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# フォークト要素のクリープコンプライアンスの計算式
def calcCreepComp(E, tau, t):
    creepComp = (1 - np.exp(-t/tau))/E
    return creepComp

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
        select = int(input('Selection (creep compliance (linear): 0, creep compliance (log): 1): '))
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
        insMod = float(input('Enter instantaneous modulus value (MPa) (default = 200.0 MPa): '))*10**6
    except ValueError:
        insMod = 200*10**6
    try:
        eta = float(input('Enter viscosity value for single Newtonian element (kPa s) (default = 1000 kPa s): '))*10**3
    except ValueError:
        eta = 10**6
    try:
        numComp = int(input('Enter the number of Voigt components (default = 1): '))
    except ValueError:
        numComp = 1
    for j in range(numComp):
        try:
            Ej = float(input('Enter modulus value of Voigt component (MPa) (default = 1.0 MPa): '))*10**6
        except ValueError:
            Ej = 10**6
        E_list.append(Ej)
        try:
            etaj = float(input('Enter viscosity value of Voigt component (kPa s) (default = 100 kPa s): '))*10**3
        except ValueError:
            etaj = 10**5
        eta_list.append(etaj)
        tauj = etaj/Ej
        tau_list.append(tauj)
    return numComp, insMod, eta, E_list, eta_list, tau_list

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
    param,_ = curve_fit(loglogFit, log_time[minFit:maxFit], log_creepComp[minFit:maxFit])
    fit_relaxMod = fittedArray(log_time, param)
    fit_result = "$J$($t$) ∝ $t^{{{0:.2f}}}$".format(param[0])
    return fit_relaxMod, fit_result

'''
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
'''

if __name__=='__main__':
    # calcul1ating creep compliance
    select, time = timeAxisChoice()
    numComp, insMod, eta, E_list, eta_list, tau_list = reqParams()
    if insMod == 0:
        infMod = 1 / np.sum(1 / np.array(E_list))   # 平衡弾性率=各フォークト要素の弾性率の逆和の逆数　式(6.28)
    else:
        infMod = 1 / (np.sum(1 / np.array(E_list)) + 1/insMod)  # 平衡弾性率=各フォークト要素の弾性率および瞬間弾性率の逆和の逆数　式(6.32)
    param_text = r'($E_i$ = {0:.1f} MPa, $E_\infty$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s, {3} components)'.format(insMod/10**6, infMod/10**6, eta/10**3, numComp)
    fitting = -1
    spectrum = -1

    # クリープコンプライアンスの計算
    if eta == 0 and insMod == 0:
        creepComp = np.zeros(len(time))                     # ゲタなし
    elif eta == 0 and insMod != 0:
        creepComp = np.ones(len(time))/insMod               # ゲタとして瞬間弾性率の逆数のnp配列, 式(6.40)    
    elif eta != 0 and insMod == 0:
        creepComp = time / eta                              # ゲタとして粘性流動項t/etaのnp配列, 式(6.44)
    else:
        creepComp = time / eta + np.ones(len(time))/insMod  # ゲタとして両方の寄与
    for j in range(numComp):
        creepComp = creepComp + calcCreepComp(E_list[j], tau_list[j], time)
    log_time = np.log10(time)
    with np.errstate(divide='ignore'):
        log_creepComp = np.log10(creepComp)
    # 各フォークト要素の緩和弾性率の計算（ただし全ての要素でゲタとして粘性流動項を加える）
    creepComp_comp = np.zeros((numComp,len(time)))
    if eta == 0 and insMod == 0:
        for j in range(numComp):
            creepComp_comp[j] = calcCreepComp(E_list[j], tau_list[j], time)
    elif eta == 0 and insMod != 0:
        for j in range(numComp):
            creepComp_comp[j] = calcCreepComp(E_list[j], tau_list[j], time) + 1/insMod
    elif eta != 0 and insMod == 0:
        for j in range(numComp):
            creepComp_comp[j] = calcCreepComp(E_list[j], tau_list[j], time) + time/eta
    else:
        for j in range(numComp):
            creepComp_comp[j] = calcCreepComp(E_list[j], tau_list[j], time) + time/eta + 1/insMod
    with np.errstate(divide='ignore'):
        log_creepComp_comp = np.log10(creepComp_comp)

    if select == 0:
        x = time
        y = log_creepComp
        x_label = r'$t$ /s'
        y_label = r'log[$J$($t$) /Pa$^{{{-1}}}$]'
        xlim = [0, np.max(x)]
        ylim = [np.min(y[1:])-0.5, np.max(y)+1.0]
        label = r'Creep compliance (linear)'
        legend_loc = 'upper left'
        c = 'r'
        savefile = './png/gen_Voigt_creepComp_linear.png'

    if select == 1:
        x = log_time
        y = log_creepComp
        y_orig = creepComp
        x_label = r'log[$t$ /s]'
        y_label = r'log[$J$($t$) /Pa$^{{{-1}}}$]'
        xlim = [np.min(x)-0.5, np.max(x)+0.5]
        ylim = [np.min(y[1:])-0.5, np.max(y)+0.5]
        label = r'Creep compliance (log)'
        legend_loc = 'upper left'
        c = 'r'
        savefile = './png/gen_Voigt_creepComp_log.png'
        try:
            fitting = int(input('Selection (curve fit: 0, no curve fit: 1): '))
        except ValueError:
            fitting = 1
        if fitting == 1:
            pass
        if fitting == 0:
            fitTimes = fitTimes()
            fit_creepComp, fit_result = curveFit(time, fitTimes)
        try:
            spectrum = int(input('Selection (with spectrum: 0, without spectrum: 1): '))
        except ValueError:
            spectrum = 1
        if spectrum == 1:
            pass
        if spectrum == 0:
            x_cropped, relaxSpectrum, log_relaxSpectrum = relaxSpectrumFunc(x,y,y_orig)

    # drawing graphs
    fig, ax = plt.subplots(figsize=(8, 6), tight_layout=True)
    ax.set_title('generalized Voigt model '+param_text)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.scatter(x, y, c=c, label=label)

    for j in range(numComp):
        ax.plot(x, log_creepComp_comp[j], c='b', linewidth=0.5)

    if fitting == 1:
        pass
    if fitting == 0:
        ax.plot(x, fit_creepComp, c='b', ls=':', label='fitted Creep compliance')
        fig.text(0.2, 0.75, fit_result)
        ax.set_xlabel(x_label)

    if spectrum == 0:
        y = [rs/insMod for rs in relaxSpectrum]
        legend_loc='upper left'
        ax1 = ax.twinx()
        ax1.scatter(x_cropped, y, c='C0', ls=':', label='Relaxation Spectrum')
        ax1.set_ylim(-0.1*np.max(y),1.5*np.max(y))
        ax1.set_ylabel('Relaxation Spectrum / Eins')
        ax1.legend(loc='upper right')
        savefile = './png/genVoigt_relaxation_modulus_log_spectrum.png'

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