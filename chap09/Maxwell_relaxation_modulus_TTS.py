# relaxation modulus of Maxwell model for TTS calculation
# 260322作成（講義：高分子レオロジーのため）
# 260617修正（さらにミスを修正）

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# preset parameters 
mu = 10**26             # entanglement density [m^(-3)]
kB = 1.38 * 10**(-23)   # Boltzmann constant [J/K] 
TA = 1727               # activation temperature [K]
V0 = 10**(-3)           # viscocity at high-temp limit [Pa s]

def calcRelaxMod(E, relaxTime, t):
    relaxMod = E*np.exp(-t/relaxTime)
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
    log_tim = np.log10(tim)
    insMod = 3*mu*kB*temp
    viscosity = V0*np.exp(TA/(temp-TV))
    relaxTime = viscosity/insMod

    param_text = '($E$($T$) = {0:.5f}*$T$ MPa, $T_{{g}}$ = {1:.0f} K)'.format(3*mu*kB/(10**6), Tg)

    relaxMod = np.zeros((len(tim), len(temp)))
    for i in range(len(temp)):
        relaxMod[:, i] = calcRelaxMod(insMod[i], relaxTime[i], tim)
    out_array = np.full_like(relaxMod, fill_value=np.nan, dtype=float)
    log_relaxMod = np.log10(relaxMod, where=(relaxMod>0), out=out_array)

    # calcul1ating shift factor
    log_aT = np.log10(np.e)*TA/(temp - TV) - 15
    aT = 10**log_aT
    log_tim_aT = np.zeros((len(tim), len(temp)))
    tim_aT = np.zeros((len(tim), len(temp)))
    for i in range(len(temp)):
        log_tim_aT[:, i] = [log_t - log_aT[i] for log_t in log_tim]
        tim_aT[:, i] = [t/aT[i] for t in tim]

    # calculating strength factor
    bT = temp/Tg

    # fitting with WLF equation

    WLF_x = temp - Tg
    WLF_y = log_aT
    WLF_fit, param = curveFit(WLF_x, WLF_y)
    fit_result = "C$_1$ = {0:.2f}, C$_2$ = {1:.2f}".format(param[0],param[1])

    # drawing graphs

    fig = plt.figure(figsize=(10,6), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title('Maxwell model '+param_text)

    try:
        with_aT = int(input('Selection of aT (without aT: 0, with aT: 1): '))
    except ValueError:
        with_aT = 0

    if with_aT == 0:
        x_label = r'log[$t$ /s]'
        ax.set_xlabel(x_label)

        try:
            select = int(input('Selection (relaxation modulus (linear): 0, relaxation modulus (log): 1): '))
        except ValueError:
            select = 0

        if select == 0:
            y_label = r'$E$ /Pa'
            ax.set_ylabel(y_label)
            savefile = './png/Maxwell_relaxation_modulus_temp_linear.png'
            for i in range(len(temp)):
                y = relaxMod[:, i]
                label_txt = r"$E$ ($T$ = {0:.0f} K)".format(temp[i])
                ax.plot(log_tim, y, c=cm.jet(0.1+float(i)/12), marker="o", lw=0.5, label=label_txt)

        if select == 1:
            y_label = r"log[$E$ /Pa]"
            ax.set_ylabel(y_label)
            savefile = './png/Maxwell_relaxation_modulus_temp_log.png'
            for i in range(len(temp)):
                y = log_relaxMod[:, i]
                label_txt = r"$E$ ($T$ = {0:.0f} K)".format(temp[i])
                ax.plot(log_tim, y, c=cm.jet(0.1+float(i)/12), marker="o", lw=0.5, label=label_txt)       
            
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0,), borderaxespad=0)
        ax.grid()
        ax.set_axisbelow(True)
        fig.savefig(savefile, dpi=300)

    if with_aT == 1:
        x_label = r'log[($t$ / $a_{{T}}$) /s]'
        ax.set_xlabel(x_label)

        try:
            select = int(input('Selection (relaxation modulus (linear): 0, relaxation modulus (log): 1): '))
        except ValueError:
            select = 0

        if select == 0:
            y_label = r"$E$ /Pa"
            ax.set_ylabel(y_label)
            try:
                with_bT = int(input('Selection of bT (without bT: 0, with bT: 1): '))
            except ValueError:
                with_bT = 0
            if with_bT == 0:
                for i in range(len(temp)):
                    x = log_tim_aT[:, i]
                    y = relaxMod[:, i]
                    label_txt = "$E$ ($T$ = {0:.0f} K)".format(temp[i])
                    ax.plot(x, y, c=cm.jet(0.1+float(i)/12), marker="o", lw=0, label=label_txt)
                    savefile = './png/Maxwell_relaxation_modulus_TTS_aT_linear.png'
            if with_bT == 1:
                for i in range(len(temp)):
                    x = log_tim_aT[:, i]
                    y = relaxMod[:, i] / bT[i]
                    label_txt = "$E$ ($T$ = {0:.0f} K)".format(temp[i])
                    ax.plot(x, y, c=cm.jet(0.1+float(i)/12), marker="o", lw=0, label=label_txt)
                    savefile = './png/Maxwell_scaled_relaxation_modulus_TTS_aTbT_linear.png'

        if select == 1:
            y_label = r"log[$E$ /Pa]"
            ax.set_ylabel(y_label)
            try:
                with_bT = int(input('Selection of bT (without bT: 0, with bT: 1): '))
            except ValueError:
                with_bT = 0
            if with_bT == 0:
                for i in range(len(temp)):
                    x = log_tim_aT[:, i]
                    y = log_relaxMod[:, i]
                    label_txt = r"$E$ ($T$ = {0:.0f} K)".format(temp[i])
                    ax.plot(x, y, c=cm.jet(0.1+float(i)/12), marker="o", lw=0, label=label_txt) 
                    savefile = './png/Maxwell_relaxation_modulus_TTS_aT_log.png'    
            if with_bT == 1:
                for i in range(len(temp)):
                    x = log_tim_aT[:, i]
                    y = log_relaxMod[:, i] - np.log10(bT[i])
                    label_txt = r"$E$ ($T$ = {0:.0f} K)".format(temp[i])
                    ax.plot(x, y, c=cm.jet(0.1+float(i)/12), marker="o", lw=0, label=label_txt)  
                    savefile = './png/Maxwell_scaled_relaxation_modulus_TTS_aTbT_log.png'  
            
        ax.grid()
        ax.set_axisbelow(True)
        fig.savefig(savefile, dpi=300)

    plt.show()