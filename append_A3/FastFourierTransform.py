# Fast Fourier Transform for various functions

import numpy as np
from scipy import fft
from scipy import signal
import matplotlib.pyplot as plt
 
# two-sided FFT
def calc_fft(data, s_rate):
    fourier = fft.fft(data)
    n = data.size
    amp = np.sqrt((fourier.real**2) + (fourier.imag**2))
    amp = amp / (n/2)
    phase = np.arctan2(fourier.imag, fourier.real)
    phase = np.degrees(phase)
    freq = fft.fftfreq(n, d=2/s_rate)
    return fourier, amp, phase, freq

if __name__=='__main__':

# parameters
    s_rate = 65536
    f = 5.0
    alpha1 = 200
    alpha2 = 30
    alpha3 = 0.05
    alpha4 = 1
    sn = 0.5

    x = 2*np.arange(0, s_rate)/s_rate - 1
    s1 = 'Selection (sin: 0, cos: 1, sawtooth: 2, Gaussian: 3, Exp: 4, Lorenzian: '
    s2 = '5, Exp+Heaviside: 6, Heaviside: 7, delta: 8, noise: 9, sin + noise: 10): '
    select_text = s1 + s2    
    try:
        select = int(input(select_text))
    except ValueError:
        select = 0
    if select == 0:
        data = np.sin(2.0*np.pi*f*x)                # sine
        legend = 'sine'
        re = 0
        im = 0
        savefile = './png/FFT_sine.png'
    if select == 1:
        data = np.cos(2.0*np.pi*f*x)                # cosine
        legend = 'cosine'
        re = 0
        im = 0
        savefile = './png/FFT_cosine.png'
    if select == 2:
        data = signal.sawtooth(2*np.pi*f*x)         # sawtooth
        legend = 'sawtooth'
        re = 0
        im = 0
        savefile = './png/FFT_sawtooth.png'
    if select == 3:
        data = np.exp(-alpha1*x**2)                 # Gaussian
        legend = 'Gaussian'
        re = 1
        im = 0
        savefile = './png/FFT_Gaussian.png'
    if select == 4:
        data = np.exp(-alpha2*np.abs(x))            # exponential
        legend = 'exponential'
        re = 1
        im = 0
        savefile = './png/FFT_exponential.png'
    if select == 5:
        data = 2*alpha3/(x**2 + alpha3**2)          # Lorentzian
        legend = 'Lorentzian'
        re = 1
        im = 0
        savefile = './png/FFT_Lorentzian.png'
    if select == 6:
        data = np.heaviside(x,0.5)*np.exp(-alpha2*x)
        legend = 'exponential (causal)'
        re = 1
        im = 1
        savefile = './png/FFT_exponential_causal.png'
    if select == 7:
        data = np.heaviside(x,0.5)                  # Heaviside
        legend = 'Heaviside'
        re = 0
        im = 1
        savefile = './png/FFT_Heaviside.png'
# 以下のデルタ関数の定義は暫定的なものであり、厳密なデルタ関数の定義ではない
    if select == 8:
        dx = 10**(-5)
        dHeaviside = np.gradient(np.heaviside(x,0.5), dx)       # Dirac delta (tentative), defined as the derivative of Heaviside function
        data = dHeaviside / np.max(dHeaviside) * 1/dx /6        # 正規化のための係数は何故か6（不明点あり）            
        legend = 'Dirac delta'
        re = 1
        im = 0
        savefile = './png/FFT_Dirac_delta.png'
    if select == 9:
        data = np.random.normal(loc=0, scale=1, size=len(x))        # noise
        legend = 'Gaussian noise'
        re = 0
        im = 0
        savefile = './png/FFT_Gaussian_noise.png'
    if select == 10:
        data = np.sin(2.0*np.pi*f*x) + sn*np.random.normal(loc=0, scale=1, size=len(x))   # sine + noise
        legend = 'sine with Gaussian noise'
        re = 0
        im = 0
        savefile = './png/FFT_sine_Gaussian_noise.png'

    fourier, amp, phase, freq = calc_fft(data, s_rate)
    real = fourier.real/ (len(data) / 2)
    imag = fourier.imag/ (len(data) / 2)
    amp_s = fft.fftshift(amp)
    real_s = fft.fftshift(real)
    imag_s = fft.fftshift(imag)
    freq_s = fft.fftshift(freq)

    pref1 = 1
    pref2 = 2*np.pi

    if select == 3:
        data_r = pref1*np.sqrt(np.pi/alpha1)*np.exp(-(pref2*freq_s)**2/(4*alpha1))
    if select == 4:
        data_r = pref1*2*alpha2/((pref2*freq_s)**2 + alpha2**2)
    if select == 5:
        data_r = pref1*2*np.pi*np.exp(-alpha3*np.abs(pref2*freq_s))
    if select == 6:
        data_r = pref1*alpha2/((pref2*freq_s)**2 + alpha2**2)
        data_i = pref1*pref2*freq_s/((pref2*freq_s)**2 + alpha2**2)
    if select == 7:
        data_i = np.sqrt(2*np.pi)*pref1*pref2*freq_s/((pref2*freq_s)**2 + alpha4**2)  # 数因子について不明点あり
    if select == 8:
        data_r = pref1*np.ones_like(freq_s)

# inverse FFT
#   ifft_time = fft.ifft(fourier)
#   ax1.plot(x, ifft_time, label='time-iFFT waveform', lw=1, color='blue')
 
    xlim = 30
    ylim = 1.2*np.max(np.maximum(amp, np.maximum(real, imag)))
    dlim = 1.5*np.max(np.abs(data))

    fig = plt.figure(figsize=(12,8), tight_layout=True)

    ax1 = fig.add_subplot(221)
    ax1.set_xlabel(r'$t$ /s')
    ax1.set_ylabel(r'$f(t)$')
    ax1.set_xlim(-1.01, 1.01)
    ax1.set_ylim(-dlim, dlim)
    ax1.plot(x, data, label=legend, lw=1, color='red')
    ax1.legend(loc='upper right')
    ax1.hlines(0, -1.01, 1.01, color='black', lw=0.5, ls='dashed')
    ax1.vlines(0, -dlim, dlim, color='black', lw=0.5, ls='dashed')

    ax2 = fig.add_subplot(222)
    ax2.set_xlabel(r'$\omega$ /s$^{-1}$')
    ax2.set_ylabel(r'$|\hat{F}(\omega)|$')
    ax2.set_xlim(-xlim, xlim)
    ax2.set_ylim(-ylim, ylim)
    ax2.plot(freq_s, amp_s, label='Amplitude', lw=1, color='blue')
    ax2.legend(loc='upper right')
    ax2.hlines(0, -xlim, xlim, color='black', lw=0.5, ls='dashed')
    ax2.vlines(0, -dlim, dlim, color='black', lw=0.5, ls='dashed')

    ax3 = fig.add_subplot(223)
    ax3.set_xlabel(r'$\omega$ /s$^{-1}$')
    ax3.set_ylabel(r'Re[${\hat{F}(\omega)}$]')
    ax3.set_xlim(-xlim, xlim)
    ax3.set_ylim(-ylim, ylim)
    ax3.plot(freq_s, real_s, label='Fourier.real', lw=1, color='blue')
    ax3.hlines(0, -xlim, xlim, color='black', lw=0.5, ls='dashed')
    ax3.vlines(0, -dlim, dlim, color='black', lw=0.5, ls='dashed')
    if re == 1:
        ax3.plot(freq_s, data_r, label='Fourier.real.theory', ls='--', lw=1, color='red')
    ax3.legend(loc='upper right')

    ax4 = fig.add_subplot(224)
    ax4.set_xlabel(r'$\omega$ /s$^{-1}$')
    ax4.set_ylabel(r'Im[${\hat{F}(\omega)}$]')
    ax4.set_xlim(-xlim, xlim)
    ax4.set_ylim(-ylim, ylim)
    ax4.plot(freq_s, imag_s, label='Fourier.imag', lw=1, color='blue')
    ax4.hlines(0, -xlim, xlim, color='black', lw=0.5, ls='dashed')
    ax4.vlines(0, -dlim, dlim, color='black', lw=0.5, ls='dashed')
    if im == 1:
        ax4.plot(freq_s, data_i, label='Fourier.imag.theory', ls='--', lw=1, color='red')
    ax4.legend(loc='upper right')
 
    fig.savefig(savefile, dpi=300)

    plt.show()