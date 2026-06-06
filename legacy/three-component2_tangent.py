import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def calc_complexMod(E, k, T, f):
    numer = 1/k + T*f*(2j/2)
    denom = 1 + T*f*(2j/2)
    comMod = E*numer/denom
    strMod = comMod.real
    losMod = comMod.imag
    return strMod, losMod

def complexMod(E, k, T, freq):
    strMod = []
    losMod = []
    for f in freq:
        s, l = calc_complexMod(E, k, T, f)
        strMod.append(s)
        losMod.append(l)
    dynamicMod = [strMod, losMod]
    return dynamicMod

def reqFreqs():
    try:
        minFreq = int(input('Enter minimum frequency in log scale (default = -2): '))
    except ValueError:
        minFreq = -2
    try:
        maxFreq = int(input('Enter maximum frequency in log scale (default = 6): '))
    except ValueError:
        maxFreq = 6
    intFreq = maxFreq - minFreq
    freqs = [minFreq, maxFreq, intFreq]
    return freqs

def reqParams():
    try:
        insMod = float(input('Enter instantaneous modulus value (GPa) (default = 1 GPa): '))*10**9
    except ValueError:
        insMod = 10**9
    try:
        infMod = float(input('Enter equilibrium modulus value (MPa) (default = 0.1 MPa): '))*10**5
    except ValueError:
        infMod = 10**5
    try:
        viscosity = float(input('Enter viscosity value (kPa s) (default = 1 MPa s): '))*10**6
    except ValueError:
        viscosity = 10**6
    modulus = insMod - infMod
    relaxTime = viscosity/modulus
    k = insMod/infMod
    return insMod, k, relaxTime

if __name__=='__main__':
    # calculating dynamic Modulus
    insMod, k, relaxTime = reqParams()
    freqs = reqFreqs()
    freq = np.logspace(freqs[0], freqs[1], freqs[2]*4+1)
    strMod = complexMod(insMod, k, relaxTime, freq)[0]
    losMod = complexMod(insMod, k, relaxTime, freq)[1]
    losTan = [l/s for (s, l) in zip(strMod, losMod)]
    log_freq = [np.log10(f) for f in freq]
#    log_strMod = [np.log10(s) for s in strMod]
#    log_losMod = [np.log10(l) for l in losMod]

    param_text = '(Ei = {0:.1f} GPa, Einf = {1:.1f} MPa, T = {2:.1f} ms)'.format(insMod/10**9, insMod/(k*10**6), relaxTime*10**3)

    # drawing graphs
    fig = plt.figure(figsize=(8,5))
    fig.suptitle('SLS model II '+param_text)
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('log(frequency /s)')
    ax1.set_ylabel('loss tangent')
    ax1.set_xlim(freqs[0],freqs[1])
    ax1.set_ylim(0, 100)
    ax1.grid()
    ax1.scatter(log_freq, losTan, c='g', label='Loss Tangent')
    ax1.legend(loc='upper right')

    fig.savefig('./png/sls_model_loss_tangent.png')

    plt.show()