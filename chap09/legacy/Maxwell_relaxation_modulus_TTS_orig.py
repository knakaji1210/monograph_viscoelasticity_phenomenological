# relaxation modulus of Maxwell model for TTS calculation
# 260322作成（講義：高分子レオロジーのため）

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# preset parameters 
mu = 10**26             # entanglement density [m^(-3)]
kB = 1.38 * 10**(-23)   # Boltzmann constant [J/K] 
TA = 1727               # activation temperature [K]
V0 = 10**(-3)           # viscocity at high-temp limit [Pa s]

def calc_relaxMod(E, relaxTime, t):
    relaxMod = E*np.exp(-t/relaxTime)
    return relaxMod

def relaxMod(E, relaxTime, tim):
    relaxMod = []
    for t in tim:
        r = calc_relaxMod(E, relaxTime, t)
        relaxMod.append(r)
    relaxMod = np.array(relaxMod)
    return relaxMod

def reqTemps():
    try:
        Tg = float(input('Enter glass-transition temperature (K) (default = 200 K): '))
    except ValueError:
        Tg = 200
    try:
        temp_max = float(input('Enter highest temperature (K) (default = 300 K): '))
    except ValueError:
        temp_max = 300
    return Tg, temp_max

def WLFFit(x, C1, C2):
    return  -C1*x / (C2 + x)

def fittedArray(x_array, param):
    fitted_array = [WLFFit(num, param[0], param[1]) for num in x_array]
    return fitted_array

def curveFit(x, y):
        param,_ = curve_fit(WLFFit, x, y)
        y_fit = fittedArray(x, param)
        return y_fit, param

if __name__=='__main__':
    # calcul1ating relaxation Modulus
    Tg, temp_max = reqTemps()
    TV = Tg - 50        # Vogel temperature [K]
    deltaT = temp_max - Tg
    temp = np.linspace(Tg, temp_max, int(deltaT/10 + 1))
    tim = np.logspace(-1, 1, 16)
    log_tim = [np.log10(t) for t in tim]
    insMod = []
    relaxTime = []
    for T in temp:
        E = 3*mu*kB*T
        V = V0*np.exp(TA/(T-TV))
        tau = V/E
        insMod.append(E)
        relaxTime.append(tau)

    param_text = '($E$($T$) = {0:.5f}*$T$ MPa, $T_{{g}}$ = {1:.0f} K)'.format(3*mu*kB/(10**6), Tg)

    relaxMod_list = []
    log_relaxMod_list = []
    for i in range(len(temp)):
        relaxationMod = relaxMod(insMod[i], relaxTime[i], tim)
        log_relaxMod = [np.log10(r) for r in relaxationMod]
        relaxMod_list.append(relaxationMod)
        log_relaxMod_list.append(log_relaxMod)  

    # calcul1ating shift factor
    log_tim_aT_list = []
    tim_aT_list = []
    log_aT = [np.log10(np.e)*TA/(T - TV) - 15 for T in temp]
    aT = [10**log for log in log_aT]
    for i in range(len(temp)):
        log_tim_aT = [log_t - log_aT[i] for log_t in log_tim]
        log_tim_aT_list.append(log_tim_aT)
        tim_aT = [t/aT[i] for t in tim]
        tim_aT_list.append(tim_aT)

    # calculating strength factor
    bT = [T/Tg for T in temp]

    # fitting with WLF equation

    wlf_x = [T - Tg for T in temp]
    wlf_y = log_aT
    wlf_fit, param = curveFit(wlf_x, wlf_y)
    fit_result = "C1 = {0:.2f}, C2 = {1:.2f}".format(param[0],param[1])

    # drawing graphs

    try:
        with_aT = int(input('Selection of aT (without aT: 0, with aT: 1): '))
    except ValueError:
        with_aT = 0

    if with_aT == 0:
        fig = plt.figure(figsize=(8,6), tight_layout=True)
        ax = fig.add_subplot(121)
        ax.set_title('Maxwell model '+param_text)

        try:
            select = int(input('Selection (relaxation modulus (linear): 0, relaxation modulus (log): 1): '))
        except ValueError:
            select = 0

        if select == 0:
            x_label = '$t$ /s'
            ax.set_xlabel(x_label)
            y_label = "$E$ /Pa"
            ax.set_ylabel(y_label)
            savefile = './png/Maxwell_relaxation_modulus_t_linear.png'
            for i in range(len(temp)):
                y = relaxMod_list[i]
                label_txt = "$E$ ($T$ = {0:.0f} K)".format(temp[i])
                ax.plot(tim, y, c=cm.jet(0.1+float(i)/12), marker="o", lw=0.5, label=label_txt)

        if select == 1:
            x_label = 'log[$t$ /s]'
            ax.set_xlabel(x_label)
            y_label = "log[$E$ /Pa]"
            ax.set_ylabel(y_label)
            savefile = './png/Maxwell_relaxation_modulus_t_log.png'
            for i in range(len(temp)):
                y = log_relaxMod_list[i]
                label_txt = "$E$ ($T$ = {0:.0f} K)".format(temp[i])
                ax.plot(log_tim, y, c=cm.jet(0.1+float(i)/12), marker="o", lw=0.5, label=label_txt)       
            
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0,), borderaxespad=0)
        ax.grid()
        ax.set_axisbelow(True)
        fig.savefig(savefile)

    if with_aT == 1:
        fig = plt.figure(figsize=(8,10), tight_layout=True)
        ax = fig.add_subplot(111)
        ax.set_title('Maxwell model '+param_text)

        try:
            select = int(input('Selection (relaxation modulus (linear): 0, relaxation modulus (log): 1): '))
        except ValueError:
            select = 0

        if select == 0:
#            x_label = '($t$ / $a_{{T}}$) /s'  # グラフが見にくいので使わない
            x_label = 'log[($t$ / $a_{{T}}$) /s]'
            ax.set_xlabel(x_label)
            y_label = "$E$ /Pa"
            ax.set_ylabel(y_label)
            try:
                with_bT = int(input('Selection of bT (without bT: 0, with bT: 1): '))
            except ValueError:
                with_bT = 0
            if with_bT == 0:
                for i in range(len(temp)):
#                    x = tim_aT_list[i]
                    x = log_tim_aT_list[i]
                    y = relaxMod_list[i]
                    label_txt = "$E$ ($T$ = {0:.0f} K)".format(temp[i])
                    ax.plot(x, y, c=cm.jet(0.1+float(i)/12), marker="o", lw=0, label=label_txt)
                    savefile = './png/Maxwell_relaxation_modulus_ts_linear.png'
            if with_bT == 1:
                for i in range(len(temp)):
#                    x = tim_aT_list[i]
                    x = log_tim_aT_list[i]
                    y = [r/bT[i] for r in relaxMod_list[i]]
                    label_txt = "$E$ ($T$ = {0:.0f} K)".format(temp[i])
                    ax.plot(x, y, c=cm.jet(0.1+float(i)/12), marker="o", lw=0, label=label_txt)
                    savefile = './png/Maxwell_scaled_relaxation_modulus_ts_linear.png'

        if select == 1:
            x_label = 'log[($t$ / $a_{{T}}$) /s]'
            ax.set_xlabel(x_label)
            y_label = "log[$E$ /Pa]"
            ax.set_ylabel(y_label)
            try:
                with_bT = int(input('Selection of bT (without bT: 0, with bT: 1): '))
            except ValueError:
                with_bT = 0
            if with_bT == 0:
                for i in range(len(temp)):
                    x = log_tim_aT_list[i]
                    y = log_relaxMod_list[i]
                    label_txt = "$E$ ($T$ = {0:.0f} K)".format(temp[i])
                    ax.plot(x, y, c=cm.jet(0.1+float(i)/12), marker="o", lw=0, label=label_txt) 
                    savefile = './png/Maxwell_relaxation_modulus_ts_log.png'    
            if with_bT == 1:
                for i in range(len(temp)):
                    x = log_tim_aT_list[i]
                    y = [log_r - np.log10(bT[i]) for log_r in log_relaxMod_list[i]]
                    label_txt = "$E$ ($T$ = {0:.0f} K)".format(temp[i])
                    ax.plot(x, y, c=cm.jet(0.1+float(i)/12), marker="o", lw=0, label=label_txt)  
                    savefile = './png/Maxwell_scaled_relaxation_modulus_ts_log.png'  
            
        ax.grid()
        ax.set_axisbelow(True)

        ax2 = fig.add_subplot(212)
        ax2.set_title('WLF equation')
        wlfx_label = '($T$ - $T_{{g}}$) /K'
        ax2.set_xlabel(wlfx_label)
        wlfy_label = 'log($a_{{T}}$)'
        ax2.set_ylabel(wlfy_label)
        ax2.scatter(wlf_x, wlf_y, c='r')
        ax2.plot(wlf_x, wlf_fit, c='b', ls='--')
        fig.text(0.6, 0.3, fit_result)

        fig.savefig(savefile, dpi=300)

    plt.show()