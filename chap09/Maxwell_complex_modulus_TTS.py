# complex modulus of Maxwell model for TTS calculation

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# preset parameters 
mu = 10**26             # entanglement density [m^(-3)]
kB = 1.38 * 10**(-23)   # Boltzmann constant [J/K] 
TA = 1727               # activation temperature [K]
V0 = 10**(-3)           # viscocity at high-temp limit [Pa s]

def calc_complexMod(E, relaxTime, f):
    numer = E*relaxTime*f*(2j/2)
    denom = 1 + relaxTime*f*(2j/2)
    comMod = numer/denom
    strMod = comMod.real
    losMod = comMod.imag
    return strMod, losMod

def complexMod(E, T, freq):
    strMod = []
    losMod = []
    for f in freq:
        s, l = calc_complexMod(E, T, f)
        strMod.append(s)
        losMod.append(l)
    dynamicMod = [strMod, losMod]
    return dynamicMod

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
    freq = np.logspace(-1, 1, 16)
    insMod = []
    relaxTime = []
    for T in temp:
        E = 3*mu*kB*T
        V = V0*np.exp(TA/(T-TV))
        tau = V/E
        insMod.append(E)
        relaxTime.append(tau)

    param_text = '(E = {0:.5f}*T MPa, Tg = {1:.0f} K)'.format(3*mu*kB/(10**6), Tg)

    log_freq = [np.log10(f) for f in freq]
    strMod_list = []
    losMod_list = []
    losTan_list = []
    log_strMod_list = []
    log_losMod_list = []
    for i in range(len(temp)):
        strMod = complexMod(insMod[i], relaxTime[i], freq)[0]
        losMod = complexMod(insMod[i], relaxTime[i], freq)[1]
        losTan = [l/s for (s, l) in zip(strMod, losMod)]
        log_strMod = [np.log10(s) for s in strMod]
        log_losMod = [np.log10(l) for l in losMod]
        strMod_list.append(strMod)
        losMod_list.append(losMod)
        losTan_list.append(losTan)
        log_strMod_list.append(log_strMod)
        log_losMod_list.append(log_losMod)    

    # calcul1ating shift factor
    log_freq_aT_list = []
    log_aT = [np.log10(np.e)*TA/(T - TV) - 15 for T in temp]
    for i in range(len(temp)):
        log_freq_aT = [log_f + log_aT[i] for log_f in log_freq]
        log_freq_aT_list.append(log_freq_aT)

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
        fig = plt.figure(figsize=(15,6), tight_layout=True)
        ax = fig.add_subplot(121)
        ax.set_title('Maxwell model '+param_text)
        freq_label = 'log[freq]'
        ax.set_xlabel(freq_label)

        try:
            select = int(input('Selection (complex modulus (linear): 0, complex modulus (log): 1, loss tangent: 2): '))
        except ValueError:
            select = 0

        if select == 0:
            y_label = "E', "+'E" /Pa'
            ax.set_ylabel(y_label)
            savefile = './png/Maxwell_complex_modulus_t_linear.png'
            for i in range(len(temp)):
                y1 = strMod_list[i]
                y2 = losMod_list[i]
                label1 = "E' (T = {0:.0f} K)".format(temp[i])
                label2 = 'E" (T = {0:.0f} K)'.format(temp[i])
                ax.plot(log_freq, y1, c=cm.jet(0.1+float(i)/12), marker="o", lw=0.5, label=label1)
                ax.plot(log_freq, y2, c=cm.jet(0.1+float(i)/12), marker="^", lw=0.5, label=label2)  

        if select == 1:
            y_label = "log[E', "+'E" /Pa]'
            ax.set_ylabel(y_label)
            savefile = './png/Maxwell_complex_modulus_t_log.png'
            for i in range(len(temp)):
                y1 = log_strMod_list[i]
                y2 = log_losMod_list[i]
                label1 = "E' (T = {0:.0f} K)".format(temp[i])
                label2 = 'E" (T = {0:.0f} K)'.format(temp[i])
                ax.plot(log_freq, y1, c=cm.jet(0.1+float(i)/12), marker="o", lw=0.5, label=label1)
                ax.plot(log_freq, y2, c=cm.jet(0.1+float(i)/12), marker="^", lw=0.5, label=label2)         
            
        if select == 2:
            y_label = 'loss tangent /'
            ax.set_ylabel(y_label)
            savefile = './png/Maxwell_loss_tangent_t.png'
            for i in range(len(temp)):
                y = losTan_list[i]
                label = "tan(d) (T = {0:.0f} K)".format(temp[i])
                ax.plot(log_freq, y, c=cm.jet(0.1+float(i)/12), marker="s", lw=0.5, label=label)
            
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0,), borderaxespad=0)
        ax.grid()
        ax.set_axisbelow(True)
        fig.savefig(savefile)

    if with_aT == 1:
        fig = plt.figure(figsize=(8,10), tight_layout=True)
        ax = fig.add_subplot(211)
        ax.set_title('Maxwell model '+param_text)
        freq_label = 'log[freq*aT]'
        ax.set_xlabel(freq_label)

        try:
            select = int(input('Selection (complex modulus (linear): 0, complex modulus (log): 1, loss tangent: 2): '))
        except ValueError:
            select = 0

        if select == 0:
            y_label = "E', "+'E" /Pa'
            ax.set_ylabel(y_label)
            savefile = './png/Maxwell_complex_modulus_ts_linear.png'
            try:
                with_bT = int(input('Selection of bT (without bT: 0, with bT: 1): '))
            except ValueError:
                with_bT = 0
            if with_bT == 0:
                for i in range(len(temp)):
                    x = log_freq_aT_list[i]
                    y1 = strMod_list[i]
                    y2 = losMod_list[i]
                    label1 = "E' (T = {0:.0f} K)".format(temp[i])
                    label2 = 'E" (T = {0:.0f} K)'.format(temp[i])
                    ax.plot(x, y1, c=cm.jet(0.1+float(i)/12), marker="o", lw=0, label=label1)
                    ax.plot(x, y2, c=cm.jet(0.1+float(i)/12), marker="^", lw=0, label=label2)  
            if with_bT == 1:
                for i in range(len(temp)):
                    x = log_freq_aT_list[i]
                    y1 = [s/bT[i] for s in strMod_list[i]]
                    y2 = [l/bT[i] for l in losMod_list[i]]
                    label1 = "E' (T = {0:.0f} K)".format(temp[i])
                    label2 = 'E" (T = {0:.0f} K)'.format(temp[i])
                    ax.plot(x, y1, c=cm.jet(0.1+float(i)/12), marker="o", lw=0, label=label1)
                    ax.plot(x, y2, c=cm.jet(0.1+float(i)/12), marker="^", lw=0, label=label2)  


        if select == 1:
            y_label = "log[E', "+'E" /Pa]'
            ax.set_ylabel(y_label)
            savefile = './png/Maxwell_complex_modulus_ts_log.png'
            try:
                with_bT = int(input('Selection of bT (without bT: 0, with bT: 1): '))
            except ValueError:
                with_bT = 0
            if with_bT == 0:
                for i in range(len(temp)):
                    x = log_freq_aT_list[i]
                    y1 = log_strMod_list[i]
                    y2 = log_losMod_list[i]
                    label1 = "E' (T = {0:.0f} K)".format(temp[i])
                    label2 = 'E" (T = {0:.0f} K)'.format(temp[i])
                    ax.plot(x, y1, c=cm.jet(0.1+float(i)/12), marker="o", lw=0, label=label1)
                    ax.plot(x, y2, c=cm.jet(0.1+float(i)/12), marker="^", lw=0, label=label2)     
            if with_bT == 1:
                for i in range(len(temp)):
                    x = log_freq_aT_list[i]
                    y1 = [log_s - np.log10(bT[i]) for log_s in log_strMod_list[i]]
                    y2 = [log_l - np.log10(bT[i]) for log_l in log_losMod_list[i]]
                    label1 = "E' (T = {0:.0f} K)".format(temp[i])
                    label2 = 'E" (T = {0:.0f} K)'.format(temp[i])
                    ax.plot(x, y1, c=cm.jet(0.1+float(i)/12), marker="o", lw=0, label=label1)
                    ax.plot(x, y2, c=cm.jet(0.1+float(i)/12), marker="^", lw=0, label=label2)       


        if select == 2:
            y_label = 'loss tangent /'
            ax.set_ylabel(y_label)
            savefile = './png/Maxwell_loss_tangent_ts.png'
            for i in range(len(temp)):
                x = log_freq_aT_list[i]
                y = losTan_list[i]
                label = "tan(d) (T = {0:.0f} K)".format(temp[i])
                ax.plot(x, y, c=cm.jet(0.1+float(i)/12), marker="o", lw=0, label=label)
            
        ax.grid()
        ax.set_axisbelow(True)

        ax2 = fig.add_subplot(212)
        ax2.set_title('WLF equation')
        wlfx_label = 'T - Tg'
        ax2.set_xlabel(wlfx_label)
        wlfy_label = 'log(aT)'
        ax2.set_ylabel(wlfy_label)
        ax2.scatter(wlf_x, wlf_y, c='r')
        ax2.plot(wlf_x, wlf_fit, c='b', ls='--')
        fig.text(0.6, 0.3, fit_result)

        fig.savefig(savefile)

    plt.show()