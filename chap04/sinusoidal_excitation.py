#　Resoponse to sinusoidal excitation

import numpy as np
import matplotlib.pyplot as plt

def reqParams():
    try:
        Freq = float(input('Enter frequency value (Hz) (default = 1.0 Hz ): '))
    except ValueError:
        Freq = 1.0
    angFreq = 2*np.pi*Freq
    try:
        select = int(input('Selection (general : 0, Hooke: 1, Newton: 2, Maxwell: 3, SLS II: 4): '))
    except ValueError:
        select = 0

    if select == 0:
        try:
            strMod = float(input('Enter storage modulus value (MPa) (default = 1.0 MPa): '))*10**6
        except ValueError:
            strMod = 10**6
        try:
            losMod = float(input('Enter loss modulus value (MPa) (default = 0.5 MPa): '))*10**6
        except ValueError:
            losMod = 0.5*10**6
        comMod = strMod + (2j/2)*losMod
        savefile = './png/sinusoidal_excitation_general.png'

    if select == 1:
        try:
            strMod = float(input('Enter modulus value (MPa) (default = 1.0 MPa): '))*10**6
        except ValueError:
            strMod = 10**6
        losMod = 0
        comMod = strMod + (2j/2)*losMod
        savefile = './png/sinusoidal_excitation_Hooke.png'    

    if select == 2:
        try:
            viscosity = float(input('Enter viscosity value (kPa s) (default = 10.0 kPa s): '))*10**3
        except ValueError:
            viscosity = 10**5
            losMod = angFreq*viscosity
        strMod = 0
        comMod = strMod + (2j/2)*losMod
        savefile = './png/sinusoidal_excitation_Newton.png'
    
    if select == 3:
        try:
            insMod = float(input('Enter instantaneous modulus value (MPa) (default = 1.0 MPa): '))*10**6
        except ValueError:
            insMod = 10**6
        try:
            viscosity = float(input('Enter viscosity value (kPa s) (default = 100.0 kPa s): '))*10**3
        except ValueError:
            viscosity = 10**5
        relaxTime = viscosity/insMod
        numer = insMod*relaxTime*angFreq*(2j/2)
        denom = 1 + relaxTime*angFreq*(2j/2)
        comMod = numer/denom
        savefile = './png/sinusoidal_excitation_Maxwell_(f={0:.1f}Hz).png'.format(Freq)

    if select == 4:
        try:
            insMod = float(input('Enter instantaneous modulus value (MPa) (default = 10.0 MPa): '))*10**6
        except ValueError:
            insMod = 10**7
        try:
            infMod = float(input('Enter equilibrium modulus value (MPa) (default = 1.0 MPa): '))*10**6
        except ValueError:
            infMod = 10**6
        try:
            viscosity = float(input('Enter viscosity value (kPa s) (default = 900.0 kPa s): '))*10**3
        except ValueError:
            viscosity = 9*10**5
        modulus = insMod - infMod
        retardTime = viscosity/modulus
        k = insMod/infMod
        numer = insMod*(1/k + retardTime*angFreq*(2j/2))
        denom = 1 + retardTime*angFreq*(2j/2)
        comMod = numer/denom
        savefile = './png/sinusoidal_excitation_slsII.png'

    return angFreq, comMod, savefile

def calc_compStrain(angFreq):
    strain_amp = 1
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
    angFreq, comMod, savefile = reqParams()
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

    fig.text(0.15, 0.15, '$f$ = {0:.3f} Hz'.format(angFreq/(2*np.pi)))
    fig.savefig(savefile, dpi=300)

    plt.show()