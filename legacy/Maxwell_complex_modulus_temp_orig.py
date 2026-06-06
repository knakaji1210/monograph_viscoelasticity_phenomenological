# complex modulus of Maxwell model

import numpy as np
import matplotlib.pyplot as plt

# preset parameters 
mu = 10**26             # entanglement density [m^(-3)]
kB = 1.38 * 10**(-23)   # Boltzmann constant [J/K] 
TA = 1727               # activation temperature [K]
V0 = 10**(-3)           # viscocity at high-temp limit [Pa s]

def calc_complexMod(E, T, f):
    numer = E*T*f*(2j/2)
    denom = 1 + T*f*(2j/2)
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
        T = V/E
        insMod.append(E)
        relaxTime.append(T)

    param_text = '(E = {0:.5f}*T MPa, Tg = {1:.0f} K)'.format(3*mu*kB/(10**6), Tg)

    log_freq = [np.log10(f) for f in freq]
    freq_label = 'log[freq]'
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

    # drawing graphs
    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(121)
    ax.set_title('Maxwell model '+param_text)
    ax.set_xlabel(freq_label)

    try:
        select = int(input('Selection (complex modulus (linear): 0, complex modulus (log): 1, loss tangent: 2): '))
    except ValueError:
        select = 0

    if select == 0:
        y_label = "E', "+'E" /Pa'
        ax.set_ylabel(y_label)
        c1 = 'r'
        c2 = 'b'
        savefile = './png/Maxwell_complex_modulus_t_linear.png'
        for i in range(len(temp)):
            y1 = strMod_list[i]
            y2 = losMod_list[i]
            label1 = "E' (T = {0:.0f} K)".format(temp[i])
            label2 = 'E" (T = {0:.0f} K)'.format(temp[i])
            ax.plot(log_freq, y1, c=c1, marker="o", lw=0.2*(i+1), label=label1)
            ax.plot(log_freq, y2, c=c2, marker="o", lw=0.2*(i+1), label=label2)  

    if select == 1:
        y_label = "log[E', "+'E" /Pa]'
        ax.set_ylabel(y_label)
        c1 = 'r'
        c2 = 'b'
        savefile = './png/Maxwell_complex_modulus_t_log.png'
        for i in range(len(temp)):
            y1 = log_strMod_list[i]
            y2 = log_losMod_list[i]
            label1 = "E' (T = {0:.0f} K)".format(temp[i])
            label2 = 'E" (T = {0:.0f} K)'.format(temp[i])
            ax.plot(log_freq, y1, c=c1, marker="o", lw=0.2*(i+1), label=label1)
            ax.plot(log_freq, y2, c=c2, marker="o", lw=0.2*(i+1), label=label2)         
            
    if select == 2:
        y_label = 'loss tangent /'
        ax.set_ylabel(y_label)
        c = 'g'
        savefile = './png/Maxwell_loss_tangent_t.png'
        for i in range(len(temp)):
            y = losTan_list[i]
            label = "tan(d) (T = {0:.0f} K)".format(temp[i])
            ax.plot(log_freq, y, c='g', marker="o", lw=0.2*(i+1), label=label)
            
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0,), borderaxespad=0)
    ax.grid()
    ax.set_axisbelow(True)
    fig.savefig(savefile)

    plt.show()