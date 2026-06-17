# complex modulus of Maxwell model for TTS calculation
# 260322修正（講義：高分子レオロジーのため）
# 260616修正（さらにミスを修正）

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# preset parameters 
mu = 10**26             # entanglement density [m^(-3)]
kB = 1.38 * 10**(-23)   # Boltzmann constant [J/K] 
TA = 1727               # activation temperature [K]
V0 = 10**(-3)           # viscocity at high-temp limit [Pa s]

def calcComplexMod(E, relaxTime, af):
    numer = E*relaxTime*af*(2j/2)
    denom = 1 + relaxTime*af*(2j/2)
    comMod = numer/denom
    strMod = comMod.real
    losMod = comMod.imag
    return strMod, losMod

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
    # calcul1ating dynamic Modulus and loss tangent
    Tg, temp_max = reqTemps()
    TV = Tg - 50        # Vogel temperature [K]
    deltaT = temp_max - Tg
    temp = np.linspace(Tg, temp_max, int(deltaT/10 + 1))
    angFreq = np.logspace(-1, 1, 16)
    insMod = 3*mu*kB*temp
    viscosity = V0*np.exp(TA/(temp-TV))
    relaxTime = viscosity/insMod

    param_text = '($E$ = {0:.5f}*$T$ MPa, $T_{{g}}$ = {1:.0f} K)'.format(3*mu*kB/(10**6), Tg)

    log_angFreq = np.log10(angFreq)
    strMod = np.zeros((len(angFreq), len(temp)))
    losMod = np.zeros((len(angFreq), len(temp)))
    for i in range(len(temp)):
        strMod[:, i] = calcComplexMod(insMod[i], relaxTime[i], angFreq)[0]
        losMod[:, i] = calcComplexMod(insMod[i], relaxTime[i], angFreq)[1]
    losTan = losMod/strMod
    log_strMod = np.log10(strMod)
    log_losMod = np.log10(losMod)  

    # calcul1ating shift factor
    log_aT = np.log10(np.e)*TA/(temp - TV) - 15
    log_afaT = np.zeros((len(angFreq), len(temp)))
    for i in range(len(temp)):
        log_afaT[:, i] = [log_af + log_aT[i] for log_af in log_angFreq]

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
        x_label = r'log($\omega\tau$)'
        ax.set_xlabel(x_label)

        try:
            select = int(input('Selection (complex modulus (linear): 0, complex modulus (log): 1, loss tangent: 2): '))
        except ValueError:
            select = 0

        if select == 0:
            y_label = r'$E^\prime$, $E^{{\prime\prime}}$ /Pa'
            ax.set_ylabel(y_label)
            savefile = './png/Maxwell_complex_modulus_temp_linear.png'
            for i in range(len(temp)):
                y1 = strMod[:, i]
                y2 = losMod[:, i]
                label1 = r'$E^\prime$ ($T$ = {0:.0f} K)'.format(temp[i])
                label2 = r'$E^{{\prime\prime}}$ ($T$ = {0:.0f} K)'.format(temp[i])
                ax.plot(log_angFreq, y1, c=cm.jet(0.1+float(i)/12), marker="o", lw=0.5, label=label1)
                ax.plot(log_angFreq, y2, c=cm.jet(0.1+float(i)/12), marker="^", lw=0.5, label=label2)  

        if select == 1:
            y_label = r'log($E^\prime$, $E^{{\prime\prime}}$ /Pa)'
            ax.set_ylabel(y_label)
            savefile = './png/Maxwell_complex_modulus_temp_log.png'
            for i in range(len(temp)):
                y1 = log_strMod[:, i]
                y2 = log_losMod[:, i]
                label1 = r'$E^\prime$ ($T$ = {0:.0f} K)'.format(temp[i])
                label2 = r'$E^{{\prime\prime}}$ ($T$ = {0:.0f} K)'.format(temp[i])
                ax.plot(log_angFreq, y1, c=cm.jet(0.1+float(i)/12), marker="o", lw=0.5, label=label1)
                ax.plot(log_angFreq, y2, c=cm.jet(0.1+float(i)/12), marker="^", lw=0.5, label=label2)         
            
        if select == 2:
            y_label = r'tan $\delta$ /'
            ax.set_ylabel(y_label)
            savefile = './png/Maxwell_loss_tangent_temp.png'
            for i in range(len(temp)):
                y = losTan[:, i]
                label = r'tan $\delta$ ($T$ = {0:.0f} K)'.format(temp[i])
                ax.plot(log_angFreq, y, c=cm.jet(0.1+float(i)/12), marker="s", lw=0.5, label=label)
            
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0,), borderaxespad=0)
        ax.grid()
        ax.set_axisbelow(True)
        fig.savefig(savefile, dpi=300)

    if with_aT == 1:
        x_label = r'log($\omega a_{{T}} \tau$)'
        ax.set_xlabel(x_label)

        try:
            select = int(input('Selection (complex modulus (linear): 0, complex modulus (log): 1, loss tangent: 2, WLF fit: 3): '))
        except ValueError:
            select = 0

        if select == 0:
            y_label = r'$E^\prime$, $E^{{\prime\prime}}$ /Pa'
            ax.set_ylabel(y_label)
            try:
                with_bT = int(input('Selection of bT (without bT: 0, with bT: 1): '))
            except ValueError:
                with_bT = 0
            if with_bT == 0:
                savefile = './png/Maxwell_complex_modulus_TTS_aT_linear.png'
                for i in range(len(temp)):
                    x = log_afaT[:, i]
                    y1 = strMod[:, i]
                    y2 = losMod[:, i]
                    label1 = r'$E^\prime$ ($T$ = {0:.0f} K)'.format(temp[i])
                    label2 = r'$E^{{\prime\prime}}$ ($T$ = {0:.0f} K)'.format(temp[i])
                    ax.plot(x, y1, c=cm.jet(0.1+float(i)/12), marker="o", lw=0, label=label1)
                    ax.plot(x, y2, c=cm.jet(0.1+float(i)/12), marker="^", lw=0, label=label2)
                    
            if with_bT == 1:
                savefile = './png/Maxwell_complex_modulus_TTS_aTbT_linear.png'
                for i in range(len(temp)):
                    x = log_afaT[:, i]
                    y1 = strMod[:, i]/bT[i]
                    y2 = losMod[:, i]/bT[i]
                    label1 = r'$E^\prime$ ($T$ = {0:.0f} K)'.format(temp[i])
                    label2 = r'$E^{{\prime\prime}}$ ($T$ = {0:.0f} K)'.format(temp[i])
                    ax.plot(x, y1, c=cm.jet(0.1+float(i)/12), marker="o", lw=0, label=label1)
                    ax.plot(x, y2, c=cm.jet(0.1+float(i)/12), marker="^", lw=0, label=label2)  
                    
        if select == 1:
            y_label = r'log($E^\prime$, $E^{{\prime\prime}}$ /Pa)'
            ax.set_ylabel(y_label)
            try:
                with_bT = int(input('Selection of bT (without bT: 0, with bT: 1): '))
            except ValueError:
                with_bT = 0
            if with_bT == 0:
                savefile = './png/Maxwell_complex_modulus_TTS_aT_log.png'
                for i in range(len(temp)):
                    x = log_afaT[:, i]
                    y1 = log_strMod[:, i]
                    y2 = log_losMod[:, i]
                    label1 = r'$E^\prime$ ($T$ = {0:.0f} K)'.format(temp[i])
                    label2 = r'$E^{{\prime\prime}}$ ($T$ = {0:.0f} K)'.format(temp[i])
                    ax.plot(x, y1, c=cm.jet(0.1+float(i)/12), marker="o", lw=0, label=label1)
                    ax.plot(x, y2, c=cm.jet(0.1+float(i)/12), marker="^", lw=0, label=label2)     
            if with_bT == 1:
                savefile = './png/Maxwell_complex_modulus_TTS_aTbT_log.png'  
                for i in range(len(temp)):
                    x = log_afaT[:, i]
                    y1 = log_strMod[:, i] - np.log10(bT[i])
                    y2 = log_losMod[:, i] - np.log10(bT[i])
                    label1 = r'$E^\prime$ ($T$ = {0:.0f} K)'.format(temp[i])
                    label2 = r'$E^{{\prime\prime}}$ ($T$ = {0:.0f} K)'.format(temp[i])
                    ax.plot(x, y1, c=cm.jet(0.1+float(i)/12), marker="o", lw=0, label=label1)
                    ax.plot(x, y2, c=cm.jet(0.1+float(i)/12), marker="^", lw=0, label=label2)
                           
        if select == 2:
            y_label = r'tan $\delta$ /'
            ax.set_ylabel(y_label)
            savefile = './png/Maxwell_loss_tangent_TTS.png'
            for i in range(len(temp)):
                x = log_afaT[:, i]
                y = losTan[:, i]
                label = r'tan $\delta$ ($T$ = {0:.0f} K)'.format(temp[i])
                ax.plot(x, y, c=cm.jet(0.1+float(i)/12), marker="o", lw=0, label=label)

        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0,), borderaxespad=0)        

        if select == 3:
            ax.set_title('WLF equation') 
            x_label = r'($T - T_{{g}}$) /K'
            ax.set_xlabel(x_label)
            y_label = r'log($a_{{T}}$)'
            ax.set_ylabel(y_label)
            label1 = r'log($a_{{T}}$)'
            label2 = r'WLF fit'
            ax.scatter(WLF_x, WLF_y, c='r', label=label1, zorder=2)
            ax.plot(WLF_x, WLF_fit, c='b', ls='--', label=label2, zorder=1)
            ax.legend(loc="upper right")
            fig.text(0.65, 0.8, fit_result)
            savefile = './png/Maxwell_WLF.png' 

        ax.grid()
        ax.set_axisbelow(True)
        fig.savefig(savefile, dpi=300)

    plt.show()