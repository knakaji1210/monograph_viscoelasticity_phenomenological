# relaxation modulus (SLS I & SLS II)

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def reqParams(model):
    # 変数の設定
    try:
        E1 = float(input('modulus 1 [MPa] (default = 1.0 MPa): '))*10**6
    except ValueError:
        E1 = 10**6                  # [Pa] 弾性率
    try:
        E2 = float(input('modulus 2 [MPa] (default = 0.2 MPa): '))*10**6
    except ValueError:
        E2 = 2*10**5                # [Pa] 弾性率
    try:
        eta = float(input('viscosity [kPa s] (default = 500.0 kPa s): '))*10**3
    except ValueError:
        eta = 5*10**5               # [Pa s] 粘度

    if model == 1:                  # SLS I
        # パラメータの計算
        insMod = E1                 # [Pa] 瞬間弾性率
        infMod = E1*E2/(E1+E2)      # [Pa] 緩和弾性率
        relaxTime = eta/E2          # [s] 緩和時間
        k = insMod/infMod
        idx_x = 1/k
        model_text = r'SLS I model '
        save_text = r'SLS1_'
        res_text = r'$E_i$ = {0:.2f} MPa, $E_\infty$ = {1:.2f} MPa, $\tau$ = {2:.2f} s, $\tau$ / $k$ = {3:.2f} s'.format(insMod/10**6, infMod/10**6, relaxTime, relaxTime/k)

    elif model == 2:                # SLS II
        # 複素弾性率の計算
        insMod = E1+E2              # [Pa] 瞬間弾性率
        infMod = E2                 # [Pa] 緩和弾性率
        relaxTime = eta/E1          # [s] 緩和時間
        k = insMod/infMod
        idx_x = 1
        model_text = r'SLS II model '
        save_text = r'SLS2_'
        res_text = r'$E_i$ = {0:.2f} MPa, $E_\infty$ = {1:.2f} MPa, $\tau$ = {2:.2f} s'.format(insMod/10**6, infMod/10**6, relaxTime)

    return E1, E2, eta, insMod, infMod, relaxTime, k, idx_x, model_text, save_text, res_text

def relaxMod(model, E, k, tau, t_elapsed):
    if model == 1:
        relaxMod = np.where(t_elapsed >= 0, E*(1/k + (1 - 1/k)*np.exp(-k*t_elapsed/tau)), 0)
        # E must be insMod
    elif model == 2:
        relaxMod = np.where(t_elapsed >= 0, E*(1/k + (1 - 1/k)*np.exp(-t_elapsed/tau)), 0)
        # E must be insMod
    return relaxMod

def relaxFunc(model, k, tau, t_elapsed):
    if model == 1:
        relaxFunc = np.where(t_elapsed >= 0, (1 - 1/k)*(1 - np.exp(-k*t_elapsed/tau)), 0)
    elif model == 2:
        relaxFunc = np.where(t_elapsed >= 0, (1 - 1/k)*(1 - np.exp(-t_elapsed/tau)), 0)
    return relaxFunc


def timeAxis(retardTime):
    log_relaxT = np.log10(retardTime)
    linearTime = np.linspace(-retardTime*0.5, retardTime*3, 400)
    scaledLinearTime = linearTime/retardTime
    logTime = np.logspace(log_relaxT-2.0, log_relaxT+2.0, 400)
    scaledLogTime = np.log10(logTime/retardTime)
    return linearTime, scaledLinearTime, logTime, scaledLogTime

def fitTimes():
    try:
        minTime = float(input('Enter minimum time for fitting in log scale (default = -2.0): '))
    except ValueError:
        minTime = -2.0
    try:
        maxTime = float(input('Enter maximum time for fitting in log scale (default = -1.5): '))
    except ValueError:
         maxTime = -1.5
    fitTimes = [minTime, maxTime]
    return fitTimes

def getNearestIdx(list, num):
    idx = np.abs(np.asarray(list) - num).argmin()
    return idx

def fitRegion(time, minNum, maxNum):
    minFit = getNearestIdx(time, minNum)
    maxFit = getNearestIdx(time, maxNum)
    fitRegion = [minFit, maxFit]
    return fitRegion

def loglogFit(x, a, b):
    return  a*x + b

def fittedArray(x_array, param):
    fitted_array = [loglogFit(num, param[0], param[1]) for num in x_array]
    return fitted_array

def curveFit(x, y, fitTimes):
        minFit = fitRegion(x, fitTimes[0], fitTimes[1])[0]
        maxFit = fitRegion(x, fitTimes[0], fitTimes[1])[1]
        param,_ = curve_fit(loglogFit, x[minFit:maxFit], y[minFit:maxFit])
        y_fit = fittedArray(x, param)
        return y_fit, param

if __name__=='__main__':
    try:
        model = int(input('Selection (SLS I : 1, SLS II: 2): '))
    except ValueError:
        model = 1    
    # calculating relaxation modulus and relaxation funcion
    E1, E2, eta, insMod, infMod, relaxTime, k, idx_x, model_text, save_text, res_text = reqParams(model)
    param_text = r'($E_1$ = {0:.1f} MPa, $E_2$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s)'.format(E1/10**6, E2/10**6, eta/10**3)
    timeAxis = timeAxis(relaxTime)
    try:
        select = int(input('Selection (relaxation modulus: 0, relaxation function: 1): '))
    except ValueError:
        select = 0

    try:
        fitting = int(input('Selection (curve fit: 0, no curve fit: 1): '))
    except ValueError:
        fitting = 0

    x1 = timeAxis[0]
    x1_scaled = timeAxis[1]
    idx1 = np.abs(x1_scaled - idx_x).argmin()
    x1_label = r'$t$/$\tau$'
    x2 = timeAxis[2]
    x2_scaled = timeAxis[3]
    idx2 = np.abs(x2_scaled - np.log(idx_x)).argmin()
    x2_label = r'log[$t$/$\tau$]'

    if select == 0:
        y1 = relaxMod(model, insMod, k, relaxTime, x1)  # E = insMod
        y1 /= 10**6             # rescale to MPa
        y1_label = r'$E$($t$) /MPa'
        label1 = r'Rleaxtion modulus (linear)'
        relaxMod = relaxMod(model, insMod, k, relaxTime, x2)
        y2 = np.log10(relaxMod)                 # creepCompの出力をnp.ndarrayに変更したのでこのように簡単に書ける
        y2_label = r'log[ $E(t)$ /Pa]'
        y_lim = [np.min(y2)-0.2, np.max(y2)+0.2]
        label2 = r'Relaxation modulus (log)'
        res_pos = [0.3, 0.8]
        if fitting == 1:
            pass
        elif fitting == 0:
            fitTimes = fitTimes()
            y_fit, param = curveFit(x2_scaled, y2, fitTimes)        
            fit_result = r'$E(t) \propto (t/\tau)^{{{0:.2f}}}$'.format(param[0])
            label_fit = r'Fitted relaxation time'
        legend_loc='upper right'
        savefile = './png/'+save_text+'relaxation modulus.png'

    elif select == 1:
        y1 = relaxFunc(model, k, relaxTime, x1)
        y1_label = r'$\phi$($t$)'
        label1 = 'Relaxation function (linear)'
        relaxFunc = relaxFunc(model, k, relaxTime, x2)
        y2 = np.log10(relaxFunc)                 # relaxFuncの出力をnp.ndarrayに変更したのでこのように簡単に書ける
        y2_label = r'log[$\phi$($t$)]'
        y_lim = [np.min(y2)-1.0, np.max(y2)+1.0]
        label2 = 'Relaxation function (log)'
        res_pos = [0.45, 0.2]
        if fitting == 1:
            pass
        elif fitting == 0:
            fitTimes = fitTimes()
            y_fit, param = curveFit(x2_scaled, y2, fitTimes)        
            fit_result = r'$\phi(t) \propto (t/\tau)^{{{0:.2f}}}$'.format(param[0])
            label_fit = r'Fitted relaxation function'
        legend_loc='upper left'
        savefile = './png/'+save_text+'relaxation_func.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,10), tight_layout=True)
    ax1 = fig.add_subplot(211)
    ax1.set_title(model_text+param_text)
    ax1.set_xlabel(x1_label)
    ax1.set_ylabel(y1_label)
    ax1.set_xlim(x1_scaled[0], x1_scaled[-1])
    ax1.set_ylim(-np.max(y1)*0.05, np.max(y1)*1.1)
    ax1.grid()
    ax1.plot(x1_scaled, y1, c='r', label=label1)
    ax1.legend(loc=legend_loc)
    ax1.vlines(idx_x, 0, y1[idx1], color='k', ls='dashed')
    ax1.set_axisbelow(True)
    ax1.text(res_pos[0], res_pos[1], res_text, transform=ax1.transAxes)
    ax2 = fig.add_subplot(212)
    ax2.set_xlabel(x2_label)
    ax2.set_ylabel(y2_label)
    ax2.set_xlim(x2_scaled[0], x2_scaled[-1])
    ax2.set_ylim(y_lim[0], y_lim[1])
    ax2.grid()
    ax2.plot(x2_scaled, y2, c='b', label=label2)
    ax2.legend(loc=legend_loc)
    ax2.vlines(np.log(idx_x), 0, y2[idx2], color='k', ls='dashed')
    ax2.set_axisbelow(True)
    if fitting == 1:
        pass
    elif fitting == 0:
        ax2.plot(x2_scaled, y_fit, c='b', ls=':', label=label_fit)
        fig.text(0.7, 0.30, fit_result)

    fig.savefig(savefile)

    plt.show()