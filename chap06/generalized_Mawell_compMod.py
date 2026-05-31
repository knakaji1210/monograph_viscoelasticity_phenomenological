# complex modulus of generalized Maxwell model (multiple freqs)
 
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# マクスウェル要素の複素弾性率の計算式
def calcComplexMod(E1, E2, tau, af):
    # E1: Maxwell spring component, E2: infMod
    numer = E1*tau*af*(2j/2)
    denom = 1 + tau*af*(2j/2)
    comMod = E2 + numer/denom
    strMod = comMod.real
    losMod = comMod.imag
    return strMod, losMod

def reqAngFreqs():
    try:
        minAngFreq = int(input('Enter minimum angular frequency in log scale (default = -3): '))
    except ValueError:
        minAngFreq = -3
    try:
        maxAngFreq = int(input('Enter maximum angular frequency in log scale (default = 6): '))
    except ValueError:
        maxAngFreq = 6
    intAngFreq = maxAngFreq - minAngFreq
    angFreqInfo = [minAngFreq, maxAngFreq, intAngFreq]
    return angFreqInfo

def angFreqAxis():
    angFreqInfo = reqAngFreqs()
    angFreq_sparse = np.logspace(angFreqInfo[0], angFreqInfo[1], angFreqInfo[2]*10+1)
    angFreq_dense = np.logspace(angFreqInfo[0], angFreqInfo[1], angFreqInfo[2]*50+1)
    angFreqAxis = [angFreq_sparse, angFreq_dense]
    return angFreqAxis

def angFreqAxisChoice():
    try:
        select = int(input('Selection (complex modulus (linear): 0, complex modulus (log): 1, loss tangent: 2, Cole-Cole plot: 3): '))
    except ValueError:
        select = 0
    if select == 0 or select == 1 or select == 2:
        angFreq = angFreqAxis()[0]
    if select == 3:
        angFreq = angFreqAxis()[1]
    return select, angFreq

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
    for i in range(numComp):
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
    return numComp, infMod, E_list, eta_list,tau_list

def fitAngFreqs():
    try:
        minAngFreq = float(input('Enter minimum angular frequency for fitting (storage) (default = 0.5): '))
    except ValueError:
        minAngFreq = 0.5
    try:
        maxAngFreq = float(input('Enter maximum angular frequency for fitting (storage) (default = 3): '))
    except ValueError:
         maxAngFreq = 3
    fitAngFreqs = [minAngFreq, maxAngFreq]
    return fitAngFreqs

def getNearestIdx(list, num):
    idx = np.abs(np.asarray(list) - num).argmin()
    return idx

def fitRegion(angFreq, minNum, maxNum):
    minFit = getNearestIdx(angFreq, 10**minNum)
    maxFit = getNearestIdx(angFreq, 10**maxNum)
    fitRegion = [minFit, maxFit]
    return fitRegion

def loglogFit(x, a, b):
    return  a*x + b

def fittedArray(x_array, param):
    fitted_array = [loglogFit(num, param[0], param[1]) for num in x_array]
    return fitted_array

def curveFit(angFreq, fitAngFreqs):
    minFit = fitRegion(angFreq, fitAngFreqs[0], fitAngFreqs[1])[0]
    maxFit = fitRegion(angFreq, fitAngFreqs[0], fitAngFreqs[1])[1]
    param,_ = curve_fit(loglogFit, log_angFreq[minFit:maxFit], log_strMod[minFit:maxFit])
    fit_strMod = fittedArray(log_angFreq, param)
    fit_result = r"$E^\prime \propto (\omega\tau)^{{{0:.2f}}}$".format(param[0])
    return fit_strMod, fit_result

if __name__=='__main__':
    # calcul1ating complex modulus and loss tangent
    select, angFreq = angFreqAxisChoice()
    numComp, infMod, E_list, eta_list,tau_list = reqParams()
    insMod = infMod + np.sum(E_list)
    param_text = r'($E_i$ = {0:.1f} MPa, $E_\infty$ = {1:.1f} MPa, {2} components)'.format(insMod/10**6, infMod/10**6, numComp)
    fitting = -1

    # 複素弾性率の計算
    strMod = infMod * np.ones(len(angFreq))      # ゲタとして平衡弾性率のnp配列を用意
    losMod = np.zeros(len(angFreq))              # ゲタとして0のnp配列を用意
    for j in range(numComp):   # 各マクスウェル要素分を追加（ただし、平衡男性率は追加済みなのでcalcComplexModの該当引数は0にする）
        strMod = strMod + calcComplexMod(E_list[j], 0, tau_list[j], angFreq)[0]
        losMod = losMod + calcComplexMod(E_list[j], 0, tau_list[j], angFreq)[1]
    scaled_strMod = strMod/insMod
    scaled_losMod = losMod/insMod
    losTan = losMod/strMod
    log_angFreq = np.log10(angFreq)
    log_strMod = np.log10(strMod)
    log_losMod = np.log10(losMod)
    # 各マクスウェル要素の緩和弾性率の計算（ただし全ての要素でゲタとして平衡弾性率を加える）
    strMod_comp = np.zeros((numComp,len(angFreq)))
    losMod_comp = np.zeros((numComp,len(angFreq)))
    for j in range(numComp):
        strMod_comp[j] = calcComplexMod(E_list[j], infMod, tau_list[j], angFreq)[0]
        losMod_comp[j] = calcComplexMod(E_list[j], infMod, tau_list[j], angFreq)[1]
    scaled_strMod_comp = strMod_comp/insMod
    scaled_losMod_comp = losMod_comp/insMod
    losTan_comp = losMod_comp/strMod_comp
    log_strMod_comp = np.log10(strMod_comp)
    log_losMod_comp = np.log10(losMod_comp)

    if select == 0:
        x = log_angFreq
        y1 = strMod/10**6
        y2 = losMod/10**6
        x_label = r'log[$\omega$]'
        y_label = r'$E^\prime$, $E^{{\prime\prime}}$ /MPa'
        xlim = [np.min(x)-0.5, np.max(x)+0.5]
        ylim = [np.max(y1)*(-0.1), np.max(y1)*1.1]
        label1 = r'Storage modulus (linear)'
        label2 = r'Loss modulus (linear)'
        c1 = 'r'
        c2 = 'b'
        a = 1
        legend_loc='upper left'
        savefile = './png/gen_Maxwell_compMod_linear.png'

    if select == 1:
        x = log_angFreq
        y1 = log_strMod
        y2 = log_losMod
        x_label = r'log[$\omega$]'
        y_label = r"log[$E^\prime$, $E^{{\prime\prime}}$ /Pa]"
        xlim = [np.min(x)-0.5, np.max(x)+0.5]
        ylim = [np.min(y1)-0.5, np.max(y1)+0.5]
        label1 = 'Storage modulus (log)'
        label2 = 'Loss modulus (log)'
        c1 = 'r'
        c2 = 'b'
        a = 1
        legend_loc='upper left'
        savefile = './png/gen_Maxwell_compMod_log.png'
        try:
            fitting = int(input('Selection (curve fit: 0, no curve fit: 1): '))
        except ValueError:
            fitting = 0
        if fitting == 1:
            pass
        if fitting == 0:
            fitAngFreqs = fitAngFreqs()
            fit_strMod, fit_result = curveFit(angFreq, fitAngFreqs)

    if select == 2:
        x = log_angFreq
        y1 = losTan
        y2 = np.zeros(len(y1))
        x_label = r'log[$\omega$]'
        y_label = r'loss tangent /'
        xlim = [np.min(x)-0.5, np.max(x)+0.5]
        ylim = [-0.1, 2.5]
        label1 = 'Loss tangent'
        label2 = ''
        c1 = 'g'
        c2 = 'b'
        a = 0
        legend_loc='upper right'
        savefile = './png/gen_Maxwell_losTan.png'

    if select == 3:
        x = scaled_strMod
        y1 = scaled_losMod
        y2 = np.zeros(len(y1))
        x_label = r"$E^\prime$ / $E_i$"
        y_label = r"$E^{{\prime\prime}}$ / $E_i$"
        xlim = [-0.05*np.max(x), 1.05*np.max(x)]
        ylim = [-0.05*np.max(y1), 1.1*np.max(y1)]
        label1 = 'Cole-Cole plot'
        label2 = ''
        c1 = 'C0'
        c2 = 'b'
        a = 0
        legend_loc='upper right'
        savefile = './png/gen_Maxwell_cole_cole.png'

    # drawing graphs
    fig, ax = plt.subplots(figsize=(8, 6), tight_layout=True)
    ax.set_title('generalized Maxwell model '+param_text)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.scatter(x, y1, c=c1, label=label1)
    ax.scatter(x, y2, c=c2, label=label2, alpha=a)

    if select == 0:
        for j in range(numComp):
            ax.plot(x, strMod_comp[j]/10**6, c='r', linewidth=0.5)
            ax.plot(x, losMod_comp[j]/10**6, c='b', linewidth=0.5)

    if select == 1:
        for j in range(numComp):
            ax.plot(x, log_strMod_comp[j], c='r', linewidth=0.5)
            ax.plot(x, log_losMod_comp[j], c='b', linewidth=0.5)

    if select == 2:
        for j in range(numComp):
            ax.plot(x, losTan_comp[j], c='g', linewidth=0.5)

    if select == 3:
        for j in range(numComp):
            ax.plot(scaled_strMod_comp[j], scaled_losMod_comp[j], c='C0', linewidth=0.5)

    if fitting == 1:
        pass
    if fitting == 0:
        ax.plot(x, fit_strMod, c='b', ls=':', label='fitted Storage modulus')
        fig.text(0.5, 0.8, fit_result)
        ax.set_xlabel(x_label)

    ax.set_xlim(xlim[0], xlim[1]) 
    ax.set_ylim(ylim[0], ylim[1])   
    ax.legend(loc=legend_loc) 
    ax.grid()
    ax.set_axisbelow(True)

    param_table = ax.table(
        cellText=[['{0:.1f}'.format(np.log10(E_list[j])) for j in range(numComp)],
                  ['{0:.1f}'.format(eta_list[j]/10**3) for j in range(numComp)],
                  ['{0:.2f}'.format(np.log10(1/tau_list[j])) for j in range(numComp)]],
        rowLabels=[r"log($E_j$ /Pa)", r"$\eta_j$ /kPa s", r"log(1/$\tau_j$ /s$^{{{-1}}}$)"],
        colLabels=['#{}'.format(j+1) for j in range(numComp)],
        loc='bottom',
        bbox=[0.0, -0.45, 1.0, 0.3]  # [x, y, 幅, 高さ] で微調整
    )
    param_table.set_fontsize(11)
    param_table.scale(1, 1.5) # セルの大きさを調整
    plt.subplots_adjust(bottom=0.2)

    fig.savefig(savefile)

    plt.show()