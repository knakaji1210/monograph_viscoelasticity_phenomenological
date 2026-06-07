# complex modulus of SLS I & SLS II models
 
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
        eta = float(input('viscosity [kPa s] (default = 100.0 kPa s): '))*10**3
    except ValueError:
        eta = 10**5               # [Pa s] 粘度

    if model == 1:                  # SLS I
        # パラメータの計算
        insMod = E1                 # [Pa] 瞬間弾性率
        infMod = E1*E2/(E1+E2)      # [Pa] 緩和弾性率
        relaxTime = eta/E2          # [s] 緩和時間
        k = insMod/infMod
        idx_x = 1/k
        model_text = r'SLS I model '
        save_text = r'SLS1_'
        res_text = r'$E_i$ = {0:.2f} MPa, $E_\infty$ = {1:.2f} MPa, $\tau$ = {2:.2f} s, $t_{{1/e}}$/$\tau$ = {3:.2f}'.format(insMod/10**6, infMod/10**6, relaxTime, idx_x)

    elif model == 2:                # SLS II
        # パラメータの計算
        insMod = E1+E2              # [Pa] 瞬間弾性率
        infMod = E2                 # [Pa] 緩和弾性率
        relaxTime = eta/E1          # [s] 緩和時間
        k = insMod/infMod
        idx_x = 1
        model_text = r'SLS II model '
        save_text = r'SLS2_'
        res_text = r'$E_i$ = {0:.2f} MPa, $E_\infty$ = {1:.2f} MPa, $\tau$ = {2:.2f} s, $t_{{1/e}}$/$\tau$ = {3:.2f}'.format(insMod/10**6, infMod/10**6, relaxTime, idx_x)

    return E1, E2, eta, insMod, infMod, relaxTime, k, idx_x, model_text, save_text, res_text


def complexMod(model, E, k, tau, af):
    if model == 1:
        numer = E*(1 + tau*af*(2j/2))
        # E must be infMod
        denom = 1 + (tau/k)*af*(2j/2)
        comMod = numer/denom
        strMod = comMod.real
        losMod = comMod.imag
    elif model == 2:
        numer = E*(1/k + tau*af*(2j/2))
        # E must be insMod
        denom = 1 + tau*af*(2j/2)
        comMod = numer/denom
        strMod = comMod.real
        losMod = comMod.imag
    return strMod, losMod

def freqAxis(relaxTime, model, k):
    if model == 1:
        centerScaledAngFreq = k
    elif model == 2:
        centerScaledAngFreq = 1
    scaledAngFreq = np.logspace(np.log10(centerScaledAngFreq)-2.0, np.log10(centerScaledAngFreq)+2.0, 51)
    angFreq = scaledAngFreq/relaxTime
    freqAxes = [angFreq, scaledAngFreq]
    return freqAxes

def fitAngFreqs():
    try:
        minAngFreq_s = float(input('Enter minimum angular frequency for fitting (storage) (default = 0.2): '))
    except ValueError:
        minAngFreq_s = 0.2
    try:
        maxAngFreq_s = float(input('Enter maximum angular frequency for fitting (storage) (default = 0.8): '))
    except ValueError:
         maxAngFreq_s = 0.8
    try:
        minAngFreq_l = float(input('Enter minimum angular frequency for fitting (loss) (default = -1.0): '))
    except ValueError:
        minAngFreq_l = -1.0
    try:
        maxAngFreq_l = float(input('Enter maximum angular frequency for fitting (loss) (default = -0.2): '))
    except ValueError:
         maxAngFreq_l = -0.2
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
    paramStr,_ = curve_fit(loglogFit, log_scaledAngFreq[minFit1:maxFit1], log_strMod[minFit1:maxFit1])
    paramLos,_ = curve_fit(loglogFit, log_scaledAngFreq[minFit2:maxFit2], log_losMod[minFit2:maxFit2])
    fit_strMod = fittedArray(log_scaledAngFreq, paramStr)
    fit_losMod = fittedArray(log_scaledAngFreq, paramLos)
    fit_result1 = r"$E^\prime \propto (\omega\tau)^{{{0:.2f}}}$".format(paramStr[0])
    fit_result2 = r"$E^{{\prime\prime}} \propto (\omega\tau)^{{{0:.2f}}}$".format(paramLos[0])
    return fit_strMod, fit_losMod, fit_result1, fit_result2

if __name__=='__main__':
    try:
        model = int(input('Selection (SLS I : 1, SLS II: 2): '))
    except ValueError:
        model = 1   
    # calcul1ating complex Modulus and loss tangent
    E1, E2, eta, insMod, infMod, relaxTime, k, idx_x, model_text, save_text, res_text = reqParams(model)
    param_text = r'($E_i$ = {0:.2f} MPa, $E_{{\infty}}$ = {1:.2f} MPa, $\tau$ = {2:.1f} s)'.format(insMod/10**6, infMod/10**6, relaxTime)
    freqAxes = freqAxis(relaxTime, model, k)
    fitting = -1
    try:
        select = int(input('Selection (complex modulus (linear): 0, complex modulus (log): 1, loss tangent: 2): '))
    except ValueError:
        select = 0

    angFreq = freqAxes[0]
    scaledAngFreq = freqAxes[1]
    scaledAngFreq_label = r'log($\omega\tau$)'
    if model == 1:
        strMod = complexMod(model, infMod, k, relaxTime, angFreq)[0]  # E = infMod
        losMod = complexMod(model, infMod, k, relaxTime, angFreq)[1]  # E = infMod
    elif model == 2:
        strMod = complexMod(model, insMod, k, relaxTime, angFreq)[0]  # E = insMod
        losMod = complexMod(model, insMod, k, relaxTime, angFreq)[1]  # E = insMod
    losTan = losMod / strMod
    log_scaledAngFreq = np.log10(scaledAngFreq)
    log_strMod = np.log10(strMod)
    log_losMod = np.log10(losMod)

    if select == 0:
        y1 = strMod/10**6             # rescale to MPa
        y2 = losMod/10**6             # rescale to MPa
        y_label = r'$E^\prime$, $E^{{\prime\prime}}$ /MPa'
        ylim = [np.max(y1)*(-0.1), np.max(y1)*1.1]
        label1 = r'$E^\prime$ (linear)'
        label2 = r'$E^{{\prime\prime}}$ (linear)'
        c1 = 'r'
        c2 = 'b'
        a = 1
        legend_loc='upper left'
        savefile = './png/'+save_text+'complex_modulus_linear.png'

    if select == 1:
        y1 = log_strMod
        y2 = log_losMod
        y_label = r'log($E^\prime$, $E^{{\prime\prime}}$ /Pa)'
        ylim = [np.min(y1)-2, np.max(y1)+1]
        label1 = r'$E^\prime$ (log)'
        label2 = r'$E^{{\prime\prime}}$ (log)'
        c1 = 'r'
        c2 = 'b'
        a = 1
        legend_loc='upper left'
        savefile = './png/'+save_text+'complex_modulus_log.png'
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
        y2 = [0 for i in range(len(y1))]
        y_label = 'tan $\delta$ /'
        ylim = [np.max(y1)*(-0.1), np.max(y1)*1.5]
        label1 = r'tan $\delta$'
        label2 = ''
        c1 = 'g'
        c2 = 'b'
        a = 0
        legend_loc='upper right'
        savefile = './png/'+save_text+'loss_tangent_compMod.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title(model_text+param_text)
    ax.set_xlabel(scaledAngFreq_label)
    ax.set_ylabel(y_label)
    ax.scatter(log_scaledAngFreq, y1, c=c1, label=label1)
    ax.scatter(log_scaledAngFreq, y2, c=c2, label=label2, alpha=a)

    if fitting == 1:
        pass
    if fitting == 0:
        ax.plot(log_scaledAngFreq, fit_strMod, c='r', ls=':', label=r'fitted $E^{\prime}$')
        ax.plot(log_scaledAngFreq, fit_losMod, c='b', ls=':', label=r'fitted $E^{{\prime\prime}}$')
        fig.text(0.2, 0.25, fit_result1)
        fig.text(0.2, 0.20, fit_result2)

    ax.set_ylim(ylim[0], ylim[1])   
    ax.legend(loc=legend_loc)
    ax.grid()
    ax.set_axisbelow(True)
    fig.savefig(savefile, dpi=300)

    plt.show()