# response of Voigt model to time-dependent stress

import numpy as np
import matplotlib.pyplot as plt

def reqParams():
    try:
        modulus = float(input('Enter modulus value (MPa) (default = 1 MPa): '))*10**6
    except ValueError:
        modulus = 10**6
    try:
        viscosity = float(input('Enter viscosity value (kPa s) (default = 10 kPa s): '))*10**3
    except ValueError:
        viscosity = 10**4
    retardTime = viscosity/modulus
    return modulus, retardTime

def func_Voigt1(modulus, retardTime):
    try:
        c = float(input('Enter c of stress = c*time (default = 0.01 MPa/ms): '))
    except ValueError:
        c = 10**7
    try:
        t1 = float(input('Enter t1 (0<=t<t1) (ms) (default = 50 ms): '))*10**(-3)
    except ValueError:
        t1 = 50*10**(-3)
    tim = np.linspace(0, t1, 400)
    stress = c*tim
    strain = (c/modulus)*(tim -retardTime * (1 - np.exp(-tim/retardTime)))
    return c, t1, tim, stress, strain

def func_Voigt2(modulus, retardTime, c, t1, e1):
    try:
        dt = float(input('Enter dt = t2 - t1 (t1<=t<t2) (default = 10 ms): '))*10**(-3)
    except ValueError:
        dt = 10*10**(-3)
    t2 = t1 + dt
    tim = np.linspace(t1, t2, 400)
    stress = c*t1*np.ones(len(tim))
    strain = (c/modulus)*(t1 -(t1 - e1/(c/modulus)) * np.exp(-(tim - t1)/retardTime))
    return t2, tim, stress, strain

def func_Voigt3(modulus, retardTime, c, t1, t2, e2):
    try:
        dt = float(input('Enter dt = t3 - t2 (t2<=t<t3) (default = 30 ms): '))*10**(-3)  
    except ValueError:
        dt = 30*10**(-3)
    t3 = t2 + dt
    tim = np.linspace(t2, t3, 400)
    stress = c*(t1 + t2 - tim)
    strain = (c/modulus)*((retardTime + t1 + t2 - tim) - (retardTime + t1 - e2/(c/modulus)) * np.exp(-(tim-t2)/retardTime))
    return t3, tim, stress, strain
    
if __name__=='__main__':
    E, T = reqParams()
    c, t1, tim1, stress1, strain1 = func_Voigt1(E, T)
    e1 = strain1[-1]    # strain at t1
    t2, tim2, stress2, strain2 = func_Voigt2(E, T, c, t1, e1)
    e2 = strain2[-1]    # strain at t2
    t3, tim3, stress3, strain3 = func_Voigt3(E, T, c, t1, t2, e2)
    tim = np.concatenate([tim1,tim2,tim3])/10**(-3)   # rescale to ms
    stress = np.concatenate([stress1,stress2,stress3])/10**6   # rescale to MPa
    strain = np.concatenate([strain1,strain2,strain3])

    param_text = r' ($E$ = {0:.1f} MPa, $\tau$ = {1:.1f} ms)'.format(E/10**6, T*10**3)

    try:
        select = int(input('Selection (strain&stress vs time: 0, stress vs strain: 1): '))
    except ValueError:
        select = 0

    if select == 0:
        fig = plt.figure(figsize=(8,10), tight_layout=True)
        ax1 = fig.add_subplot(211)
        ax1.set_title('Voigt model for time-dependent strain'+param_text)
        ax1.set_xlabel(r'$t$ /ms')
        ax1.set_ylabel(r'$\sigma$ /MPa')
        ax1.set_xlim(0, t3*10**3)
        ax1.set_ylim(0, np.max(stress)*1.2)
        ax1.grid()
        ax1.set_axisbelow(True)
        ax1.plot(tim, stress, c='r', lw=2, label='Time-dependent stress')
        ax1.legend(loc='upper right')

        ax2 = fig.add_subplot(212)
        ax2.set_xlabel(r'$t$ /ms')
        ax2.set_ylabel(r'$\epsilon$ /')
        ax2.set_xlim(0, t3*10**3)
        ax2.set_ylim(0, np.max(strain)*1.2)
        ax2.grid()
        ax2.set_axisbelow(True)
        ax2.plot(tim, strain, c='b', lw=2, label='Response to time-dependent stress')
        ax2.legend(loc='upper right')
        savefile = './png/Voigt_time-dependent_strain_(tau={0:.1f}ms).png'.format(T*10**3)

    elif select == 1:
        fig = plt.figure(figsize=(8,5), tight_layout=True)
        ax1 = fig.add_subplot(111)
        ax1.set_title('Voigt model for time-dependent strain'+param_text)
        ax1.set_xlabel(r'$\epsilon$ /')
        ax1.set_ylabel(r'$\sigma$ /MPa')
        ax1.set_xlim(0, np.max(strain)*1.2)
        ax1.set_ylim(0, np.max(stress)*1.2)
        ax1.grid()
        ax1.set_axisbelow(True)
        ax1.plot(strain, stress, c='r', lw=2, label='stress-strain curve')
        ax1.legend(loc='upper right')
        savefile = './png/Voigt_stress-strain_curve_(tau={0:.1f}ms).png'.format(T*10**3)

    fig.savefig(savefile, dpi=300)
    plt.show()