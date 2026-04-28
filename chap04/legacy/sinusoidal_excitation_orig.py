#　Resoponse to sinusoidal excitation

import numpy as np
import matplotlib.pyplot as plt

def reqParams():
    try:
        Freq = float(input('Enter frequency value (/s) (default = 1 /s): '))
    except ValueError:
        Freq = 1
    angFreq = 2*np.pi*Freq
    try:
        select = int(input('Selection (general : 0, Hooke: 1, Newton: 2, Maxwell: 3, SLS II: 4): '))
    except ValueError:
        select = 0

    if select == 0:
        try:
            strMod = float(input('Enter storage modulus value (MPa) (default = 1 MPa): '))*10**6
        except ValueError:
            strMod = 10**6
        try:
            losMod = float(input('Enter loss modulus value (MPa) (default = 0.5 MPa): '))*10**6
        except ValueError:
            losMod = 0.5*10**6
        comMod = strMod + (2j/2)*losMod

    if select == 1:
        try:
            strMod = float(input('Enter modulus value (MPa) (default = 1 MPa): '))*10**6
        except ValueError:
            strMod = 10**6
        losMod = 0
        comMod = strMod + (2j/2)*losMod    

    if select == 2:
        try:
            viscosity = float(input('Enter viscosity value (kPa s) (default = 10 kPa s): '))*10**3
        except ValueError:
            viscosity = 10**5
            losMod = angFreq*viscosity
        strMod = 0
        comMod = strMod + (2j/2)*losMod
    
    if select == 3:
        try:
            insMod = float(input('Enter instantaneous modulus value (MPa) (default = 1 MPa): '))*10**6
        except ValueError:
            insMod = 10**6
        try:
            viscosity = float(input('Enter viscosity value (kPa s) (default = 100 kPa s): '))*10**3
        except ValueError:
            viscosity = 10**5
        relaxTime = viscosity/insMod
        numer = insMod*relaxTime*angFreq*(2j/2)
        denom = 1 + relaxTime*angFreq*(2j/2)
        comMod = numer/denom

    if select == 4:
        try:
            insMod = float(input('Enter instantaneous modulus value (MPa) (default = 10 MPa): '))*10**6
        except ValueError:
            insMod = 10**7
        try:
            infMod = float(input('Enter equilibrium modulus value (MPa) (default = 1 MPa): '))*10**6
        except ValueError:
            infMod = 10**6
        try:
            viscosity = float(input('Enter viscosity value (kPa s) (default = 900 kPa s): '))*10**3
        except ValueError:
            viscosity = 9 * 10**5
        modulus = insMod - infMod
        retardTime = viscosity/modulus
        k = insMod/infMod
        numer = insMod*(1/k + retardTime*angFreq*(2j/2))
        denom = 1 + retardTime*angFreq*(2j/2)
        comMod = numer/denom

    return angFreq, comMod

def calc_compStrain(angFreq):
    strain_amp = 1
    time_min = 0
    period = 2*np.pi/angFreq
    time_max = period*3
    time = np.linspace(time_min, time_max, 200)
    timeInfo = [time_min, time_max, time]
    omegaTime = [angFreq*t for t in time]
    compStrain = [strain_amp*np.exp((2j/2)*wt) for wt in omegaTime]

    return timeInfo, compStrain

def calc_compStress(compStrain, comMod):
    compStress = [comMod * e for e in compStrain]
    
    return compStress


if __name__=='__main__':
    angFreq, comMod = reqParams()
    timeInfo, compStrain = calc_compStrain(angFreq)
    time_min = timeInfo[0]
    time_max = timeInfo[1]
    time = timeInfo[2]
    strain = [e.real for e in compStrain]
    compStress = calc_compStress(compStrain, comMod)
    stress = [s.real for s in compStress]
    norm_stress = [s/abs(comMod) for s in stress]
    param_text = "(E' = {0:.2f} MPa, ".format(comMod.real/10**6) + 'E" = {0:.2f} MPa)'.format(comMod.imag/10**6)
    
    fig = plt.figure(tight_layout=True)
    ax1 = fig.add_subplot(111)
    ax1.set_title('Response to sinusoidal excitation '+param_text)
    ax1.set_xlim(time_min, time_max)
    ax1.set_ylim(1.5*np.min(strain),1.5*np.max(strain))
    ax1.set_xlabel('t /s')
    ax1.set_ylabel('Strain')
    ax1.plot(time, strain, c='b', label='Strain')

    ax2 = ax1.twinx()
    ax2.plot(time, norm_stress, c='r', label='Stress/Modulus')
    ax2.set_ylim(2.0*np.min(norm_stress),2.0*np.max(norm_stress))
    ax2.set_ylabel('Stress/Modulus')

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()

    ax1.legend(h1+h2, l1+l2, loc='upper right')
    ax1.grid()

    fig.text(0.15, 0.15, 'freq = {0:.3f} Hz'.format(angFreq/(2*np.pi)))
    fig.savefig('./png/sinusoidal_excitation.png')

    plt.show()