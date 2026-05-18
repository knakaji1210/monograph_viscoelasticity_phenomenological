# creep compliance (SLS I & SLS II)

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
        retardTime = eta/E2          # [s] 緩和時間
        k = insMod/infMod
        idx_x = 1
        model_text = r'SLS I model '
        save_text = r'SLS1_'

    elif model == 2:                # SLS II
        # パラメータの計算
        insMod = E1+E2              # [Pa] 瞬間弾性率
        infMod = E2                 # [Pa] 緩和弾性率
        retardTime = eta/E1          # [s] 緩和時間
        k = insMod/infMod
        idx_x = k
        model_text = r'SLS II model '
        save_text = r'SLS2_'

    return E1, E2, eta, insMod, infMod, retardTime, k, idx_x, model_text, save_text

def creepComp(model, E, k, tau, t_elapsed):
    if model == 1:
        creepComp = np.where(t_elapsed >= 0, (1 - (1 - 1/k)*np.exp(-t_elapsed/tau))/E, 0)
        # E must be infMod
    elif model == 2:
        creepComp = np.where(t_elapsed >= 0, (1 - (1 - 1/k)*np.exp(-t_elapsed/(k*tau)))/E, 0)
        # E must be infMod
    return creepComp

def creepFunc(model, k, tau, t_elapsed):
    if model == 1:
        creepFunc = np.where(t_elapsed >= 0, -1 + k*(1 - (1 - 1/k)*np.exp(-t_elapsed/tau)), 0)
    if model == 2:
        creepFunc = np.where(t_elapsed >= 0, (k - 1)*(1 - np.exp(-t_elapsed/(k*tau))), 0)
    return creepFunc

def timeAxis(retardTime, model, k):
    log_retardT = np.log10(retardTime)
    if model == 1:
        linearTime = np.linspace(-retardTime*0.5, retardTime*3, 400)
        logTime = np.logspace(log_retardT-2.0, log_retardT+2.0, 400)
    if model == 2:
        linearTime = np.linspace(-retardTime*0.5, retardTime*3*k, 400)
        logTime = np.logspace(log_retardT-2.0, log_retardT+3.0, 400)
    scaledLinearTime = linearTime/retardTime
    scaledLogTime = np.log10(logTime/retardTime)
    return linearTime, scaledLinearTime, logTime, scaledLogTime

def fitTimes():
    try:
        minTime = float(input('Enter minimum time for fitting in log scale (default = -2.0): '))
    except ValueError:
        minTime = -2.0
    try:
        maxTime = float(input('Enter maximum time for fitting in log scale (default = -1.0): '))
    except ValueError:
         maxTime = -1.0
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
    # calculating creep compliance and creep funcion
    E1, E2, eta, insMod, infMod, retardTime, k, idx_x, model_text, save_text = reqParams(model)
    param_text = r'($E_1$ = {0:.1f} MPa, $E_2$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s)'.format(E1/10**6, E2/10**6, eta/10**3)
    res_text = r'$E_i$ = {0:.2f} MPa, $E_\infty$ = {1:.2f} MPa, $\tau$ = {2:.2f} s, $t_{{1/e}}$/$\tau$ = {3:.2f}'.format(insMod/10**6, infMod/10**6, retardTime, idx_x)
    timeAxis = timeAxis(retardTime, model, k)
    try:
        select = int(input('Selection (creep compliance: 0, creep function: 1): '))
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
    idx2 = np.abs(x2_scaled - np.log10(idx_x)).argmin()
    x2_label = r'log[$t$/$\tau$]'

    if select == 0:
        y1 = creepComp(model, infMod, k, retardTime, x1)
        y1 *= 10**6             # rescale to MPa
        y1_label = r'$J$($t$) /MPa$^{{-1}}$'
        label1 = r'Creep compliance (linear)'
        creepComp = creepComp(model, infMod, k, retardTime, x2)
        y2 = np.log10(creepComp)                 # creepCompの出力をnp.ndarrayに変更したのでこのように簡単に書ける
        y2_label = r'log[ $J$($t$) /Pa$^{{-1}}$]'
        y_lim = [np.min(y2)-0.2, np.max(y2)+0.2]
        label2 = r'Creep compliance (log)'
        res_pos = [0.05, 0.85]
        if fitting == 1:
            pass
        if fitting == 0:
            fitTimes = fitTimes()
            y_fit, param = curveFit(x2_scaled, y2, fitTimes)        
            fit_result = r'$J(t) \propto (t/\tau)^{{{0:.2f}}}$'.format(param[0])
            label_fit = r'Fitted creep compliance'
        legend_loc='upper left'
        savefile = './png/'+save_text+'creep_compliance.png'

    if select == 1:
        y1 = creepFunc(model, k, retardTime, x1)
        y1_label = r'$\psi$($t$)'
        label1 = 'Creep function (linear)'
        creepFunc = creepFunc(model, k, retardTime, x2)
        y2 = np.log10(creepFunc)                 # creepFuncの出力をnp.ndarrayに変更したのでこのように簡単に書ける
        y2_label = r'log[$\psi$($t$)]'
        y_lim = [np.min(y2)-1.0, np.max(y2)+1.0]
        label2 = 'Creep function (log)'
        res_pos = [0.05, 0.85]
        if fitting == 1:
            pass
        if fitting == 0:
            fitTimes = fitTimes()
            y_fit, param = curveFit(x2_scaled, y2, fitTimes)        
            fit_result = r'$\psi(t) \propto (t/\tau)^{{{0:.2f}}}$'.format(param[0])
            label_fit = r'Fitted creep function'
        legend_loc='upper left'
        savefile = './png/'+save_text+'creep_func.png'

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
    ax2.vlines(np.log10(idx_x), 0, y2[idx2], color='k', ls='dashed')
    ax2.set_axisbelow(True)
    if fitting == 1:
        pass
    if fitting == 0:
        ax2.plot(x2_scaled, y_fit, c='b', ls=':', label=label_fit)
        ax2.text(0.05, 0.5, fit_result, transform=ax2.transAxes)
    
    fig.savefig(savefile, dpi=300)

    plt.show()