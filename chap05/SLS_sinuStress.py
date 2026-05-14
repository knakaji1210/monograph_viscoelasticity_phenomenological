#　Resoponse to sinusoidal stress (SLS I & SLS II)

import numpy as np
import matplotlib.pyplot as plt

def reqParams():
    try:
        Freq = float(input('Enter frequency value (Hz) (default = 0.3183 Hz ): '))
    except ValueError:
        Freq = 0.3183
    angFreq = 2*np.pi*Freq
    try:
        stress_amp = float(input('stress amplitute [MPa] (default = 0.04 MPa): '))*10**6
    except ValueError:
        stress_amp = 0.04*10**6     # [Pa] 応力振幅
    try:
        select = int(input('Selection (SLS I : 0, SLS II: 1): '))
    except ValueError:
        select = 0

    if select == 0:
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
        # 複素コンプライアンスの計算
        insMod = E1                 # [Pa] 瞬間弾性率
        infMod = E1*E2/(E1+E2)      # [Pa] 緩和弾性率
        k = insMod/infMod
        tau = eta/E2                # [s] 緩和時間
        numer = k + tau*angFreq*(2j/2)        
        denom = insMod*(1 + tau*angFreq*(2j/2)) 
        comComp = numer/denom
        modeltext = r'SLS I model'
        figtext = r'$\tau$ = {0:.2f} s, $f$ = {1:.3f} Hz, $\omega\tau$ = {2:.2f}'.format(tau, Freq, angFreq*tau)
        savefile = './png/sinuStress_SLS_I.png'

    if select == 1:
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
        # 複素コンプライアンスの計算
        insMod = E1+E2              # [Pa] 瞬間弾性率
        infMod = E2                 # [Pa] 緩和弾性率
        k = insMod/infMod
        tau = eta/E1                # [s] 緩和時間    
        numer = 1 + tau*angFreq*(2j/2)
        denom = insMod*(1/k + tau*angFreq*(2j/2))
        comComp = numer/denom
        modeltext = r'SLS II model'
        figtext = r'$\tau$ = {0:.2f} s, $f$ = {1:.3f} Hz, $\omega\tau$ = {2:.2f}'.format(tau, Freq, angFreq*tau)
        savefile = './png/sinuStress_SLS_II.png'

    return angFreq, stress_amp, comComp, modeltext, figtext, savefile

def calc_compStress(stress_amp, angFreq):
    time_min = 0
    period = 2*np.pi/angFreq
    time_max = period*3
    time = np.linspace(time_min, time_max, 200)
    timeInfo = [time_min, time_max, time]
    compStress = stress_amp*np.exp((2j/2)*angFreq*time)

    return timeInfo, compStress

def calc_compStrain(compStress, comComp):
    compStrain = compStress*comComp
    
    return compStrain


if __name__=='__main__':
    angFreq, stress_amp, comComp, modeltext, figtext, savefile = reqParams()
    timeInfo, compStress = calc_compStress(stress_amp, angFreq)
    time_min = timeInfo[0]
    time_max = timeInfo[1]
    time = timeInfo[2]
    stress = np.real(compStress)
    compStrain = calc_compStrain(compStress, comComp)
    strain = np.real(compStrain)
    stress /= 10**6             # rescale to MPa    
    param_text = r'($J^\prime$ = {0:.2f} MPa$^{{-1}}$, $J^{{\prime\prime}}$ = {1:.2f} MPa$^{{-1}}$)'.format(comComp.real*10**6,-comComp.imag*10**6)
    
    fig = plt.figure(tight_layout=True)
    ax1 = fig.add_subplot(111)
    ax1.set_title('Response to sinusoidal stress '+ param_text)
    ax1.set_xlim(time_min, time_max)
    ax1.set_ylim(1.5*np.min(stress),1.5*np.max(stress))
    ax1.set_xlabel(r'$t$ /s')
    ax1.set_ylabel(r'$\sigma$ /MPa')
    ax1.plot(time, stress, c='r', label=r'$\sigma$ /MPa (input)')

    ax2 = ax1.twinx()
    ax2.plot(time, strain, c='b', label=r'$\epsilon$ (output)')
    ax2.set_ylim(2.0*np.min(strain),2.0*np.max(strain))
    ax2.set_ylabel(r'$\epsilon$ /')

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()

    ax1.legend(h1+h2, l1+l2, loc='upper right')
    ax1.grid()

    fig.text(0.15, 0.85, modeltext)
    fig.text(0.15, 0.15, figtext)

    fig.savefig(savefile, dpi=300)

    plt.show()