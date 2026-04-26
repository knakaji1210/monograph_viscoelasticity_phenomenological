# relaxation modulus of Maxwell model

import numpy as np
# from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def relaxMod(E, tau, t_elapsed):
    relaxMod = np.where(t_elapsed >= 0, E * np.exp(-t_elapsed/tau), 0)
    return relaxMod

def relaxFunc(t_elapsed, tau):
    relaxFunc = np.where(t_elapsed >= 0, 1 - np.exp(-t_elapsed/tau), 0)
    return relaxFunc

def reqParams():
    try:
        insMod = float(input('Enter instantaneous modulus value (MPa) (default = 1 MPa): '))*10**6
    except ValueError:
        insMod = 10**6
    try:
        viscosity = float(input('Enter viscosity value (kPa s) (default = 100 kPa s): '))*10**3
    except ValueError:
        viscosity = 10**5
    relaxTime = viscosity/insMod
    return insMod, relaxTime

def timeAxis(relaxTime):
    log_relaxT = np.log10(relaxTime)
    linearTime = np.linspace(-relaxTime*0.5, relaxTime*3, 400)
#    scaledLinearTime = [t/relaxTime for t in linearTime]
    scaledLinearTime = linearTime/relaxTime
    logTime = np.logspace(log_relaxT-1.0, log_relaxT+1.0, 400)
    scaledLogTime = np.log10(logTime/relaxTime)
    return linearTime, scaledLinearTime, logTime, scaledLogTime

if __name__=='__main__':
    # calculating relaxation modulus and relaxation funcion
    insMod, relaxTime = reqParams()
    param_text = r'($E_{{i}}$ = {0:.1f} MPa, $\tau$ = {1:.1f} ms)'.format(insMod/10**6, relaxTime*10**3)
    timeAxis = timeAxis(relaxTime)
    try:
        select = int(input('Selection (relaxation modulus: 0, relaxation function: 1): '))
    except ValueError:
        select = 0

    x1 = timeAxis[0]
    x1_scaled = timeAxis[1]
    idx1 = np.abs(x1_scaled - 1).argmin()
    x1_label = r'$t$/$\tau$'
    x2 = timeAxis[2]
    x2_scaled = timeAxis[3]
    idx2 = np.abs(x2_scaled - 0).argmin()
    x2_label = r'log[$t$/$\tau$]'

    if select == 0:        
        y1 = relaxMod(insMod, relaxTime, x1)
        y1 /= 10**6             # rescale to MPa
        y1_label = r'$E$($t$) /MPa'
        label1 = 'Relaxation modulus (linear)'
        relaxMod = relaxMod(insMod, relaxTime, x2)
        y2 = np.log10(relaxMod)                 # relaxModの出力をnp.ndarrayに変更したのでこのように簡単に書ける
        y2_label = r'log[$E$($t$) /Pa]'
        y_lim = [np.min(y2)-1, np.max(y2)+1]
        label2 = 'Relaxation modulus (log)'
        legend_loc='upper right'
        savefile = './png/Maxwell_relax_modulus.png'

    if select == 1:
        x1 = timeAxis[0]
        x1_scaled = timeAxis[1]
        x1_label = r'$t$/$\tau$'
        y1 = relaxFunc(x1, relaxTime)
        y1_label = r'$\varphi$($t$)'
        label1 = 'Relaxation function (linear)'
        relaxFunc = relaxFunc(x2, relaxTime)
        y2 = np.log10(relaxFunc)                 # relaxFuncの出力をnp.ndarrayに変更したのでこのように簡単に書ける
        y2_label = r'log[$\varphi$($t$)]'
        y_lim = [np.min(y2)-0.2, 0.2]
        label2 = 'Relaxation function (log)'
        legend_loc='upper left'
        savefile = './png/Maxwell_relax_func.png'

    # drawing graphs
    fig = plt.figure(figsize=(8,10), tight_layout=True)
    ax1 = fig.add_subplot(211)
    ax1.set_title('Maxwell model '+param_text)
    ax1.set_xlabel(x1_label)
    ax1.set_ylabel(y1_label)
    ax1.set_xlim(x1_scaled[0], x1_scaled[-1])
    ax1.set_ylim(-np.max(y1)*0.05, np.max(y1)*1.1)
    ax1.grid()
    ax1.plot(x1_scaled, y1, c='r', label=label1)
    ax1.legend(loc=legend_loc)
    ax1.vlines(1, 0, y1[idx1], color='k', ls='dashed')
    ax1.set_axisbelow(True)
    ax2 = fig.add_subplot(212)
    ax2.set_xlabel(x2_label)
    ax2.set_ylabel(y2_label)
    ax2.set_xlim(x2_scaled[0], x2_scaled[-1])
    ax2.set_ylim(y_lim[0], y_lim[1])
    ax2.grid()
    ax2.plot(x2_scaled, y2, c='b', label=label2)
    ax2.legend(loc=legend_loc)
    ax2.vlines(0, 0, y2[idx2], color='k', ls='dashed')
    ax2.set_axisbelow(True)

    fig.savefig(savefile, dpi=300)

    plt.show()