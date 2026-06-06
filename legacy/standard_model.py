import numpy as np
import matplotlib.pyplot as plt

def calc_std_model(k0, k, c, f):
    strMod = k0 + k * c**2 * f**2 / (k**2 + c**2 * f**2)
    losMod = k**2 * c * f / (k**2 + c**2 * f**2)
    return strMod, losMod

def std_model(k0, k, c, freq):
    strMod = []
    losMod = []
    for f in freq:
        s, l = calc_std_model(k0, k, c, f)
        strMod.append(s)
        losMod.append(l)
    log_strMod = [np.log10(s) for s in strMod]
    log_losMod = [np.log10(l) for l in losMod]
    return log_strMod, log_losMod


if __name__=='__main__':
    freq = np.logspace(-4, 6, 37)
    log_freq = [np.log10(f) for f in freq]
    log_strMod1, log_losMod1 = std_model(1, 10, 100, freq)
    log_strMod2, log_losMod2 = std_model(1, 10**2, 50, freq)
    log_strMod3, log_losMod3 = std_model(1, 10**3, 30, freq)
    log_strMod4, log_losMod4 = std_model(1, 10**4, 20, freq)

    fig = plt.figure()
    plt.title('Standard Linear Viscoelastic Model')
    plt.xlabel('Frequency')
    plt.ylabel('Modulus')
    ax1 = plt.plot(log_freq,log_strMod1, c='r', label='Storage Modulus #1')
    ax2 = plt.plot(log_freq,log_losMod1, c='b', label='Loss Modulus #1')
    ax3 = plt.plot(log_freq,log_strMod2, c='r', label='Storage Modulus #2')
    ax4 = plt.plot(log_freq,log_losMod2, c='b', label='Loss Modulus #2')
    ax5 = plt.plot(log_freq,log_strMod3, c='r', label='Storage Modulus #3')
    ax6 = plt.plot(log_freq,log_losMod3, c='b', label='Loss Modulus #3')
    ax7 = plt.plot(log_freq,log_strMod4, c='r', label='Storage Modulus #4')
    ax8 = plt.plot(log_freq,log_losMod4, c='b', label='Loss Modulus #4')
    plt.legend(loc='upper left')
    plt.xlim(-4,6)
    plt.ylim(-2,6)
    plt.show()