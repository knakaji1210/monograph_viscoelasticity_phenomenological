# complex compliance of generalized Voigt model (multiple freqs)
 
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# フォークト要素の複素コンプライアンスの計算式
def calcComplexComp(E, tau, af):
    numer = E**(-1)
    denom = 1 + af*tau*(2j/2)
    comComp = numer/denom
    strComp = comComp.real
    losComp = -comComp.imag
    return strComp, losComp

def reqAngFreqs():
    try:
        minAngFreq = int(input('Enter minimum angular frequency in log scale (default = -3): '))
    except ValueError:
        minAngFreq = -3
    try:
        maxAngFreq = int(input('Enter maximum angular frequency in log scale (default = 7): '))
    except ValueError:
        maxAngFreq = 7
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
        insMod = float(input('Enter instantaneous modulus value (MPa) (default = 200.0 MPa): '))*10**6
    except ValueError:
        insMod = 200*10**6
    try:
        eta = float(input('Enter viscosity value for single Newtonian element (MPa s) (default = 1000 MPa s): '))*10**6
    except ValueError:
        eta = 10**9
    try:
        numComp = int(input('Enter the number of Voigt components (default = 1): '))
    except ValueError:
        numComp = 1
    for i in range(numComp):
        try:
            Ej = float(input('Enter modulus value of Voigt component (MPa) (default = 1 MPa): '))*10**6
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
    return numComp, insMod, eta, E_list, eta_list,tau_list

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
    param,_ = curve_fit(loglogFit, log_angFreq[minFit:maxFit], log_strComp[minFit:maxFit])
    fit_strComp = fittedArray(log_angFreq, param)
    fit_result = r"$J^\prime \propto (\omega\tau)^{{{0:.2f}}}$".format(param[0])
    return fit_strComp, fit_result

if __name__=='__main__':
    # calcul1ating complex compliance and loss tangent
    select, angFreq = angFreqAxisChoice()
    numComp, insMod, eta, E_list, eta_list, tau_list = reqParams()
    if insMod == 0:
        infMod = 1 / np.sum(1 / np.array(E_list))   # 平衡弾性率=各フォークト要素の弾性率の逆和の逆数　式(6.28)
    else:
        infMod = 1 / (np.sum(1 / np.array(E_list)) + 1/insMod)  # 平衡弾性率=各フォークト要素の弾性率および瞬間弾性率の逆和の逆数　式(6.32)
    param_text = r'($E_i$ = {0:.1f} MPa, $E_\infty$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s, {3} components)'.format(insMod/10**6, infMod/10**6, eta/10**3, numComp)
    fitting = -1

    # 複素コンプライアンスの計算
    strComp = np.zeros(len(angFreq))
    if eta == 0 and insMod == 0:
        strComp = np.zeros(len(angFreq))                # ゲタなし
        losComp = np.zeros(len(angFreq))                # ゲタなし
    elif eta == 0 and insMod != 0:
        strComp = np.ones(len(angFreq))/insMod          # ゲタとして瞬間弾性率の逆数のnp配列, 式(6.47)
        losComp = np.zeros(len(angFreq))                # ゲタなし
    elif eta != 0 and insMod == 0:
        strComp = np.zeros(len(angFreq))                # ゲタなし
        losComp = 1 / (eta * angFreq)                   # ゲタとして粘性流動項1/(ηω)のnp配列, 式(6.5X)
    else:
        strComp = np.ones(len(angFreq))/insMod          # ゲタとして瞬間弾性率の逆数のnp配列, 式(6.47)
        losComp = 1 / (eta * angFreq)                   # ゲタとして粘性流動項1/(ηω)のnp配列, 式(6.5X)
    for j in range(numComp):
        strComp = strComp + calcComplexComp(E_list[j], tau_list[j], angFreq)[0]
        losComp = losComp + calcComplexComp(E_list[j], tau_list[j], angFreq)[1]
    log_angFreq = np.log10(angFreq)
    with np.errstate(divide='ignore'):
        log_strComp = np.log10(strComp)
        log_losComp = np.log10(losComp)
    scaled_strComp = strComp*infMod
    scaled_losComp = losComp*infMod
    losTan = losComp/strComp
    log_angFreq = np.log10(angFreq)
    log_strComp = np.log10(strComp)
    log_losComp = np.log10(losComp)
    # 各フォークト要素の複素コンプライアンスの計算（ただし全ての要素でゲタとして粘性流動項を加える）
    strComp_comp = np.zeros((numComp,len(angFreq)))
    losComp_comp = np.zeros((numComp,len(angFreq)))    
    if eta == 0 and insMod == 0:
        for j in range(numComp):
            strComp_comp[j] = calcComplexComp(E_list[j], tau_list[j], angFreq)[0]                           # ゲタなし
            losComp_comp[j] = calcComplexComp(E_list[j], tau_list[j], angFreq)[1]                           # ゲタなし 
    elif eta == 0 and insMod != 0:
        for j in range(numComp):
            strComp_comp[j] = calcComplexComp(E_list[j], tau_list[j], angFreq)[0] + 1 / insMod              # 瞬間コンプライアンス項
            losComp_comp[j] = calcComplexComp(E_list[j], tau_list[j], angFreq)[1]                           # ゲタなし
    elif eta != 0 and insMod == 0:
        for j in range(numComp):
            strComp_comp[j] = calcComplexComp(E_list[j], tau_list[j], angFreq)[0]                           # ゲタなし
            losComp_comp[j] = calcComplexComp(E_list[j], tau_list[j], angFreq)[1] + 1 / (eta * angFreq)     # 粘性流動項
    else:
        for j in range(numComp):
            strComp_comp[j] = calcComplexComp(E_list[j], tau_list[j], angFreq)[0] + 1 / insMod              # 瞬間コンプライアンス項
            losComp_comp[j] = calcComplexComp(E_list[j], tau_list[j], angFreq)[1] + 1 / (eta * angFreq)     # 粘性流動項
    scaled_strComp_comp = strComp_comp*infMod
    scaled_losComp_comp = losComp_comp*infMod
    losTan_comp = losComp_comp/strComp_comp
    log_strComp_comp = np.log10(strComp_comp)
    log_losComp_comp = np.log10(losComp_comp)

    if select == 0:
        x = log_angFreq
        y1 = strComp*10**6
        y2 = losComp*10**6
        x_label = r'log[$\omega$]'
        y_label = r'$J^\prime$, $J^{{\prime\prime}}$ /MPa$^{{{-1}}}$'
        xlim = [np.min(x)-0.5, np.max(x)+0.5]
        ylim = [np.max(y1)*(-0.1), np.max(y1)*1.1]
        label1 = r'Storage compliance (linear)'
        label2 = r'Loss compliance (linear)'
        c1 = 'r'
        c2 = 'b'
        a = 1
        legend_loc='upper right'
        savefile = './png/gen_Voigt_compComp_linear.png'

    if select == 1:
        x = log_angFreq
        y1 = log_strComp
        y2 = log_losComp
        x_label = r'log[$\omega$]'
        y_label = r"log[$J^\prime$, $J^{{\prime\prime}}$ /Pa$^{{{-1}}}$]"
        xlim = [np.min(x)-0.5, np.max(x)+0.5]
        ylim = [np.min(y1)-0.5, np.max(y1)+0.5]
        label1 = 'Storage compliance (log)'
        label2 = 'Loss compliance (log)'
        c1 = 'r'
        c2 = 'b'
        a = 1
        legend_loc='upper right'
        savefile = './png/gen_Voigt_compComp_log.png'
        try:
            fitting = int(input('Selection (curve fit: 0, no curve fit: 1): '))
        except ValueError:
            fitting = 0
        if fitting == 1:
            pass
        if fitting == 0:
            fitAngFreqs = fitAngFreqs()
            fit_strComp, fit_result = curveFit(angFreq, fitAngFreqs)

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
        savefile = './png/gen_Voigt_losTan.png'

    if select == 3:
        x = scaled_strComp
        y1 = scaled_losComp
        y2 = np.zeros(len(y1))
        x_label = r"$J^\prime$ / $J_i$"
        y_label = r"$J^{{\prime\prime}}$ / $J_i$"
        xlim = [-0.05*np.max(x), 1.05*np.max(x)]
        ylim = [-0.05*np.max(y1), 1.1*np.max(y1)]
        label1 = 'Cole-Cole plot'
        label2 = ''
        c1 = 'C0'
        c2 = 'b'
        a = 0
        legend_loc='upper right'
        savefile = './png/gen_Voigt_cole_cole.png'

    # drawing graphs
    fig, ax = plt.subplots(figsize=(8, 6), tight_layout=True)
    ax.set_title('generalized Voigt model '+param_text)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.scatter(x, y1, c=c1, label=label1)
    ax.scatter(x, y2, c=c2, label=label2, alpha=a)

    if select == 0:
        for j in range(numComp):
            ax.plot(x, strComp_comp[j]*10**6, c='r', linewidth=0.5)
            ax.plot(x, losComp_comp[j]*10**6, c='b', linewidth=0.5)

    if select == 1:
        for j in range(numComp):
            ax.plot(x, log_strComp_comp[j], c='r', linewidth=0.5)
            ax.plot(x, log_losComp_comp[j], c='b', linewidth=0.5)

    if select == 2:
        for j in range(numComp):
            ax.plot(x, losTan_comp[j], c='g', linewidth=0.5)

    if select == 3:
        for j in range(numComp):
            ax.plot(scaled_strComp_comp[j], scaled_losComp_comp[j], c='C0', linewidth=0.5)

    if fitting == 1:
        pass
    if fitting == 0:
        ax.plot(x, fit_strComp, c='b', ls=':', label='fitted Storage compliance')
        fig.text(0.25, 0.75, fit_result)
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