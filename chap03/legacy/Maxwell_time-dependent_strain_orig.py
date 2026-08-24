# response of Maxwell model to time-dependent strain

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
    relaxTime = viscosity/modulus
    return modulus, relaxTime

def func_Maxwell1(modulus, relaxTime):
    try:
        c = float(input('Enter c of strain = c*time (default = 10): '))
    except ValueError:
        c = 10
    try:
        t1 = float(input('Enter t1 (0<=t<t1) (ms) (default = 50 ms): '))*10**(-3)
    except ValueError:
        t1 = 50*10**(-3)
    tim = np.linspace(0, t1, 500)
    strain = np.array([c*t for t in tim])
    stress = np.array([c*modulus*relaxTime*(1 - np.exp(-t/relaxTime)) for t in tim])
    stress = np.array(stress)
    stress /= 10**6             # rescale to MPa
    return c, t1, tim, strain, stress

def func_Maxwell2(modulus, relaxTime, c, t1):
    try:
        dt = float(input('Enter dt = t2 - t1 (t1<=t<t2) (default = 10 ms): '))*10**(-3)
    except ValueError:
        dt = 10*10**(-3)
    t2 = t1 + dt
    tim = np.linspace(t1, t2, 500)
    strain = np.array([c*t1 for t in tim])
    stress = np.array([c*modulus*relaxTime*(np.exp(t1/relaxTime) - 1)*np.exp(-t/relaxTime) for t in tim])
    stress = np.array(stress)
    stress /= 10**6             # rescale to MPa
    return t2, tim, strain, stress

def func_Maxwell3(modulus, relaxTime, c, t1, t2):
    try:
        dt = float(input('Enter dt = t3 - t2 (t2<=t<t3) (default = 30 ms): '))*10**(-3)  
    except ValueError:
        dt = 30*10**(-3)
    t3 = t2 + dt
    tim = np.linspace(t2, t3, 500)
    strain = np.array([c*(t1 + t2 -t) for t in tim])
    stress = np.array([c*modulus*relaxTime*((np.exp(t1/relaxTime) + np.exp(t2/relaxTime) - 1)*np.exp(-t/relaxTime) - 1) for t in tim])
    stress = np.array(stress)
    stress /= 10**6             # rescale to MPa
    return t3, tim, strain, stress

if __name__=='__main__':
    E, T = reqParams()
    c, t1, tim1, strain1, stress1 = func_Maxwell1(E, T)
    t2, tim2, strain2, stress2 = func_Maxwell2(E, T, c, t1)
    t3, tim3, strain3, stress3 = func_Maxwell3(E, T, c, t1, t2)
    strain = np.concatenate([strain1,strain2,strain3])
    stress = np.concatenate([stress1,stress2,stress3])
    tim1 *=10**3    # rescale to ms
    tim2 *=10**3    # rescale to ms
    tim3 *=10**3    # rescale to ms
    min_stress = np.min(stress)
    max_stress = np.max(stress)

    param_text = r' ($E$ = {0:.1f} MPa, $\tau$ = {1:.1f} ms)'.format(E/10**6, T*10**3)

    fig = plt.figure(figsize=(8,10), tight_layout=True)
    ax1 = fig.add_subplot(211)
    ax1.set_title('Maxwell model for time-dependent strain'+param_text)
    ax1.set_xlabel(r'$t$ /ms')
    ax1.set_ylabel(r'$\epsilon$ /')
    ax1.set_xlim(0, t3*10**3)
    ax1.set_ylim(np.max(strain)*(-0.2), np.max(strain)*1.2)
    ax1.grid()
    ax1.set_axisbelow(True)
    ax1.plot(tim1, strain1, c='b', lw=2, label='Time-dependent strain')
    ax1.plot(tim2, strain2, c='b', lw=2)
    ax1.plot(tim3, strain3, c='b', lw=2)
    ax1.legend(loc='upper right')

    ax2 = fig.add_subplot(212)
    ax2.set_xlabel(r'$t$ /ms')
    ax2.set_ylabel(r'$\sigma$ /MPa')
    ax2.set_xlim(0, t3*10**3)
    ax2.set_ylim(abs(min_stress)*(-1.2), abs(max_stress)*1.2)
    ax2.grid()
    ax2.set_axisbelow(True)
    ax2.plot(tim1, stress1, c='r', lw=2, label='Response to time-dependent strain')
    ax2.plot(tim2, stress2, c='r', lw=2)
    ax2.plot(tim3, stress3, c='r', lw=2)
    ax2.legend(loc='upper right')

    savefile = './png/Maxwell_time-dependent_strain_(tau={0:.1f}ms).png'.format(T*10**3)
    fig.savefig(savefile, dpi=300)

    plt.show()