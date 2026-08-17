# dynamic modulus of fractional Maxwell model
 
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os

def reqParams():
    try:
        insMod = float(input('Enter instantaneous modulus value (MPa) (default = 1 MPa): '))*10**6
    except ValueError:
        insMod = 10**6
    try:
        modulus = float(input('Enter modulus value of spring-pot (MPa) (default = 0.1 MPa): '))*10**6
    except ValueError:
        modulus = 10**5
    try:
        viscosity = float(input('Enter viscosity value of spring-pot (kPa s) (default = 100 kPa s): '))*10**3
    except ValueError:
        viscosity = 10**5
    try:
        nu = float(input('Enter fractional power (0 < nu < 1) (default = 0.5): '))
    except ValueError:
        nu = 0.5
    kappa = modulus / insMod
    retardTime = kappa**(1/nu)*viscosity/modulus
    return insMod, modulus, retardTime, nu

def freqAxes(relaxTime):
    centerAngFreq = 1 / relaxTime
    angFreq = np.logspace(int(np.log10(centerAngFreq))-5, int(np.log10(centerAngFreq))+5, 51)
    scaledAngFreq = angFreq*relaxTime
    freqAxes = [angFreq, scaledAngFreq]
    return freqAxes

# fractional_Maxwell_complex_modulus.pyとの比較して、
# この部分だけが異なる。E'とE"について、顕な式を用いて計算するように変更した。
# fractional_Maxwell_complex_modulus.pyの結果と同じになることを確認済み。
def calc_dynamicMod(E, tau, af, nu):
    x = (tau*af)**nu
    strNumer = E * (np.cos(np.pi*nu/2)*x + x**2)
    losNumer = E * np.sin(np.pi*nu/2)*x
    # E must be insMod
    denom = 1 + 2*np.cos(np.pi*nu/2)*x + x**2
    strMod = strNumer/denom
    losMod = losNumer/denom
    return strMod, losMod

def dynamicMod(E, tau, angFreq, nu):
    strMod, losMod = calc_dynamicMod(E, tau, angFreq, nu)
    dynamicMod = [strMod, losMod]
    return dynamicMod

def fitAngFreqs():
    try:
        minAngFreq_s = float(input('Enter minimum frequency for fitting (storage) (default = -5): '))
    except ValueError:
        minAngFreq_s = -5
    try:
        maxAngFreq_s = float(input('Enter maximum frequency for fitting (storage) (default = -3): '))
    except ValueError:
         maxAngFreq_s = -3
    try:
        minAngFreq_l = float(input('Enter minimum frequency for fitting (loss) (default = -5): '))
    except ValueError:
        minAngFreq_l = -5
    try:
        maxAngFreq_l = float(input('Enter maximum frequency for fitting (loss) (default = -3): '))
    except ValueError:
         maxAngFreq_l = -3
    fitFreqs = [minAngFreq_s, maxAngFreq_s, minAngFreq_l, maxAngFreq_l]
    return fitFreqs

def getNearestIdx(list, num):
    idx = np.abs(np.asarray(list) - num).argmin()
    return idx

def fitRegion(angFreq, minNum, maxNum):
    idx1 = getNearestIdx(angFreq, 10**minNum)
    idx2 = getNearestIdx(angFreq, 10**maxNum)
    # インデックスの大小関係を保証し、終端を含めるために +1 する
    minFit = min(idx1, idx2)
    maxFit = max(idx1, idx2) + 1
    fitRegion = [minFit, maxFit]
    return fitRegion

def loglogFit(x, a, b):
    return  a*x + b

def fittedArray(x_array, param):
    fitted_array = [loglogFit(num, param[0], param[1]) for num in x_array]
    return fitted_array

def curveFit(scaledAngFreq, fitAngFreqs):
    minFit1, maxFit1 = fitRegion(scaledAngFreq, fitAngFreqs[0], fitAngFreqs[1])
    minFit2, maxFit2 = fitRegion(scaledAngFreq, fitAngFreqs[2], fitAngFreqs[3])
    
    # 1. 貯蔵弾性率（E'）のフィッティング
    # 選択した範囲のデータがすべて NaN でないかチェック
    if not np.isnan(log_strMod[minFit1:maxFit1]).all() and (maxFit1 - minFit1) > 1:
        paramStr,_ = curve_fit(loglogFit, log_scaledAngFreq[minFit1:maxFit1], log_strMod[minFit1:maxFit1])
        fit_strMod = fittedArray(log_scaledAngFreq, paramStr)
        fit_result1 = r"$E^\prime \propto (\omega\tau)^{{{0:.2f}}}$".format(paramStr[0])
    else:
        # すべて NaN、またはデータ不足の場合はフィッティングをスキップ
        fit_strMod = np.full_like(log_scaledAngFreq, np.nan)
        fit_result1 = r"$E^\prime$: Fit Skipped (NaN)"

    # 2. 損失弾性率（E''）のフィッティング
    # 選択した範囲のデータがすべて NaN でないかチェック (nu=0 の時はここが True になります)
    if not np.isnan(log_losMod[minFit2:maxFit2]).all() and (maxFit2 - minFit2) > 1:
        paramLos,_ = curve_fit(loglogFit, log_scaledAngFreq[minFit2:maxFit2], log_losMod[minFit2:maxFit2])
        fit_losMod = fittedArray(log_scaledAngFreq, paramLos)
        fit_result2 = r"$E^{{\prime\prime}} \propto (\omega\tau)^{{{0:.2f}}}$".format(paramLos[0])
    else:
        # すべて NaN、またはデータ不足の場合はフィッティングをスキップ
        fit_losMod = np.full_like(log_scaledAngFreq, np.nan)
        fit_result2 = r"$E^{{\prime\prime}}$: Fit Skipped (NaN)"
    
    return fit_strMod, fit_losMod, fit_result1, fit_result2

if __name__=='__main__':
    # 保存先ディレクトリの作成（存在しない場合のエラー防止）
    os.makedirs('./png', exist_ok=True)

    # calcul1ating complex modulus and loss tangent
    insMod, modulus, relaxTime, nu = reqParams()
    param_text = r'($E_i$ = {0:.1f} MPa, $E$ = {1:.1f} MPa, $\tau$ = {2:.1f} ms, $\nu$ = {3:.2f})'.format(insMod/10**6, modulus/10**6, relaxTime*10**3, nu)
    freqAxes = freqAxes(relaxTime)
    fitting = -1
    try:
        select = int(input('Selection (complex modulus (linear): 0, complex modulus (log): 1, loss tangent: 2): '))
    except ValueError:
        select = 0

    angFreq = freqAxes[0]
    scaledAngFreq = freqAxes[1]
    scaledAngFreq_label = r'log($\omega\tau$)'
    strMod = dynamicMod(modulus, relaxTime, angFreq, nu)[0]
    losMod = dynamicMod(modulus, relaxTime, angFreq, nu)[1]
    losTan = np.divide(losMod, strMod, where=strMod!=0, out=np.full_like(strMod, np.nan))
    log_scaledAngFreq = np.log10(scaledAngFreq)
    log_strMod = np.log10(strMod, where=strMod>0, out=np.full_like(strMod, np.nan))
    log_losMod = np.log10(losMod, where=losMod>0, out=np.full_like(losMod, np.nan))

    if select == 0:
        y1 = strMod/10**6             # rescale to MPa
        y2 = losMod/10**6             # rescale to MPa
        y_label = r'$E^\prime$, $E^{{\prime\prime}}$ /MPa'
        ylim = [np.nanmax(np.concatenate([y1, y2]))*(-0.2), np.nanmax(np.concatenate([y1, y2]))*1.4]
        label1 = r'$E^\prime$ (linear)'
        label2 = r'$E^{{\prime\prime}}$ (linear)'
        c1 = 'r'
        c2 = 'b'
        a = 1
        legend_loc='upper left'
        savefile = './png/fractional_Maxwell_dynamic_modulus_linear_nu{:.2f}.png'.format(nu)

    if select == 1:
        y1 = log_strMod
        y2 = log_losMod
        y_label = r'log($E^\prime$, $E^{{\prime\prime}}$ /Pa)'
        ylim = [np.nanmin(np.concatenate([y1, y2]))-2, np.nanmax(np.concatenate([y1, y2]))+2]
        label1 = r'$E^\prime$ (log)'
        label2 = r'$E^{{\prime\prime}}$ (log)'
        c1 = 'r'
        c2 = 'b'
        a = 1
        legend_loc='upper left'
        savefile = './png/fractional_Maxwell_dynamic_modulus_log_nu{:.2f}.png'.format(nu)

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
        y2 = np.zeros(len(y1))
        y_label = 'tan $\delta$ /'
        ylim = [-0.2, 1.4]
        label1 = r'tan $\delta$'
        label2 = ''
        c1 = 'g'
        c2 = 'b'
        a = 0
        legend_loc='upper right'
        savefile = './png/fractional_Maxwell_loss_tangent_nu{:.2f}_2.png'.format(nu)

    # drawing graphs
    fig = plt.figure(figsize=(8,5), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title('fractional Maxwell model '+param_text)
    ax.set_xlabel(scaledAngFreq_label)
    ax.set_ylabel(y_label)
    ax.scatter(log_scaledAngFreq, y1, c=c1, label=label1)
    ax.scatter(log_scaledAngFreq, y2, c=c2, label=label2, alpha=a)

    # fittingを実行した（fitting == 0）ときのみプロットする
    if fitting == 0:
        ax.plot(log_scaledAngFreq, fit_strMod, c='r', lw=1, ls=':', label=r'fitted $E^{\prime}$')
        ax.plot(log_scaledAngFreq, fit_losMod, c='b', lw=1, ls=':', label=r'fitted $E^{{\prime\prime}}$')
#        fig.text(0.7, 0.40, fit_result1)
#        fig.text(0.7, 0.35, fit_result2)
        fig.text(0.7, 0.30, fit_result1)
        fig.text(0.7, 0.25, fit_result2)

    ax.set_ylim(ylim[0], ylim[1])   
    ax.legend(loc=legend_loc)
    ax.grid()
    ax.set_axisbelow(True)
    fig.savefig(savefile, dpi=300)

    plt.show()