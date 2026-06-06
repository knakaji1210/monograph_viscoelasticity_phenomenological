import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

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

def reqFreqs():
    try:
        minFreq = int(input('Enter minimum frequency in log scale (default = -1): '))
    except ValueError:
        minFreq = -1
    try:
        maxFreq = int(input('Enter maximum frequency in log scale (default = 5): '))
    except ValueError:
        maxFreq = 5
    intFreq = maxFreq - minFreq
    freqs = [minFreq, maxFreq, intFreq]
    return freqs

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

if __name__=='__main__':
    # calculating dynamic Modulus
    E, T = reqParams()
    freqs = reqFreqs()
    freq = np.logspace(freqs[0], freqs[1], freqs[2]*4+1)
    strMod = complexMod(E, T, freq)[0]
    losMod = complexMod(E, T, freq)[1]
    losTan = [l/s for (s, l) in zip(strMod, losMod)]
    log_freq = [np.log10(f) for f in freq]
#    log_strMod = [np.log10(s) for s in strMod]
#    log_losMod = [np.log10(l) for l in losMod]


    param_text = '(E = {0:.1f} MPa, T = {1:.1f} ms)'.format(E/10**6, T*10**3)

    # drawing graphs
    fig = plt.figure(figsize=(8,5))
    fig.suptitle('Maxwell model '+param_text)
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('log(frequency /s)')
    ax1.set_ylabel('loss tangent')
    ax1.set_xlim(freqs[0],freqs[1])
#    ax1.set_ylim(np.max(log_losMod)*(-0.2), np.max(log_losMod)*1.2)
    ax1.grid()
    ax1.scatter(log_freq, losTan, c='g', label='Loss Tangent')
    ax1.legend(loc='upper right')

    fig.savefig('./png/maxwell_loss_tangent.png')

    plt.show()