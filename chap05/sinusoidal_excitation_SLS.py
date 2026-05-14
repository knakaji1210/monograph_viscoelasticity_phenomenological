#　Resoponse to sinusoidal excitation (SLS I & SLS II)

import numpy as np
import matplotlib.pyplot as plt

def reqParams():
    try:
        Freq = float(input('Enter frequency value (Hz) (default = 1.0 Hz ): '))
    except ValueError:
        Freq = 1.0
    angFreq = 2*np.pi*Freq
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
            eta = float(input('viscosity [kPa s] (default = 500.0 kPa s): '))*10**3
        except ValueError:
            eta = 5*10**5               # [Pa s] 粘度
        # 複素弾性率の計算
        insMod = E1                 # [Pa] 瞬間弾性率
        infMod = E1*E2/(E1+E2)      # [Pa] 緩和弾性率
        k = insMod/infMod
        tau = eta/E2                # [s] 緩和時間        
        numer = insMod*(1 + tau*angFreq*(2j/2))
        denom = k + tau*angFreq*(2j/2)
        comMod = numer/denom
        figtext = r'$\tau$ = {0:.2f} s, $f$ = {1:.3f} Hz, $\omega\tau$ = {2:.2f}'.format(tau, Freq, angFreq*tau)
        savefile = './png/sinusoidal_excitation_SLS_I.png'

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
            eta = float(input('viscosity [kPa s] (default = 500.0 kPa s): '))*10**3
        except ValueError:
            eta = 5*10**5               # [Pa s] 粘度
        # 複素弾性率の計算
        insMod = E1+E2              # [Pa] 瞬間弾性率
        infMod = E2                 # [Pa] 緩和弾性率
        k = insMod/infMod
        tau = eta/E1                # [s] 緩和時間    
        numer = insMod*(1/k + tau*angFreq*(2j/2))
        denom = 1 + tau*angFreq*(2j/2)
        comMod = numer/denom
        figtext = r'$\tau$ = {0:.2f} s, $f$ = {1:.3f} Hz, $\omega\tau$ = {2:.2f}'.format(tau, Freq, angFreq*tau)
        savefile = './png/sinusoidal_excitation_SLS_II.png'

    return angFreq, comMod, figtext, savefile

def calc_compStrain(angFreq):
    strain_amp = 0.25
    time_min = 0
    period = 2*np.pi/angFreq
    time_max = period*3
    time = np.linspace(time_min, time_max, 200)
    timeInfo = [time_min, time_max, time]
    omegaTime = np.array([angFreq*t for t in time])
    compStrain = np.array([strain_amp*np.exp((2j/2)*wt) for wt in omegaTime])

    return timeInfo, compStrain

def calc_compStress(compStrain, comMod):
    compStress = np.array([comMod * e for e in compStrain])
    
    return compStress


if __name__=='__main__':
    angFreq, comMod, figtext, savefile = reqParams()
    timeInfo, compStrain = calc_compStrain(angFreq)
    time_min = timeInfo[0]
    time_max = timeInfo[1]
    time = timeInfo[2]
    strain = np.array([e.real for e in compStrain])
    compStress = calc_compStress(compStrain, comMod)
    stress = np.array([s.real for s in compStress])
    stress /= 10**6             # rescale to MPa    
    param_text = r'($E^\prime$ = {0:.2f} MPa, $E^{{\prime\prime}}$ = {1:.2f} MPa)'.format(comMod.real/10**6,comMod.imag/10**6)
    
    fig = plt.figure(tight_layout=True)
    ax1 = fig.add_subplot(111)
    ax1.set_title('Response to sinusoidal excitation '+param_text)
    ax1.set_xlim(time_min, time_max)
    ax1.set_ylim(1.5*np.min(strain),1.5*np.max(strain))
    ax1.set_xlabel(r'$t$ /s')
    ax1.set_ylabel(r'$\epsilon$ /')
    ax1.plot(time, strain, c='b', label=r'$\epsilon$')

    ax2 = ax1.twinx()
    ax2.plot(time, stress, c='r', label=r'$\sigma$ /MPa')
    ax2.set_ylim(2.0*np.min(stress),2.0*np.max(stress))
    ax2.set_ylabel(r'$\sigma$ /MPa')

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()

    ax1.legend(h1+h2, l1+l2, loc='upper right')
    ax1.grid()

    fig.text(0.15, 0.15, figtext)

    fig.savefig(savefile, dpi=300)

    plt.show()