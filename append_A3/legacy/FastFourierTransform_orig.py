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
    freq = fft.fftfreq(n, d=1/s_rate)
    return fourier, amp, phase, freq

if __name__=='__main__':
    s_rate = 65536
    f = 20
    x = np.arange(0, s_rate)/s_rate
#   data = np.sin(2.0*np.pi*f*x)                                # sine wave
#   data = np.cos(2.0*np.pi*f*x)                                # cosine wave
#   data = np.sin(2.0*np.pi*f*x) + 0.2*np.random.randn(len(x))  # sine wave + noise
#   data = signal.sawtooth(2*np.pi*f*x)                         # sawtooth wave
#   data = np.random.normal(loc=0, scale=1, size=len(x))
#   data = 2*np.exp(-100*x**2)
#   data = 2*np.exp(-100*(x-0.5)**2)
#   data = np.exp(-10*np.abs(x-0.5))
    data = 1/((x-0.5)**2 + 0.01)

    fourier, amp, phase, freq = calc_fft(data, s_rate)
    real = fourier.real/ (len(data) / 2)
    imag = fourier.imag/ (len(data) / 2)
    amp_s = fft.fftshift(amp)
    real_s = fft.fftshift(real)
    imag_s = fft.fftshift(imag)
    freq_s = fft.fftshift(freq)
    ifft_time = fft.ifft(fourier)
 
    fig = plt.figure(figsize=(12,8), tight_layout=True)

    ax1 = fig.add_subplot(221)
    ax1.set_xlabel('time /s')
    ax1.set_ylabel('amplitude')
    ax1.plot(x, data, label='original waveform', lw=1, color='red')
#   ax1.plot(x, ifft_time, label='time-iFFT waveform', lw=1, color='blue')

    ax2 = fig.add_subplot(222)
    ax2.set_xlabel('frequency /Hz')
    ax2.set_ylabel('Fourier (amplitude)')
    ax2.set_xlim(-100, 100)
    ax2.set_ylim(-1.2, 1.2)
    ax2.plot(freq_s, amp_s, label='Amplitude', lw=1, color='blue')

    ax3 = fig.add_subplot(223)
    ax3.set_xlabel('frequency /Hz')
    ax3.set_ylabel('Fourier (real)')
    ax3.set_xlim(-100, 100)
    ax3.set_ylim(-1.2, 1.2)
    ax3.plot(freq_s, real_s, label='Fourier.real', lw=1, color='blue')

    ax4 = fig.add_subplot(224)
    ax4.set_xlabel('frequency /Hz')
    ax4.set_ylabel('Fourier (imaginary)')
    ax4.set_xlim(-100, 100)
    ax4.set_ylim(-1.2, 1.2)
    ax4.plot(freq_s, imag_s, label='Fourier.imag', lw=1, color='blue')
 
    plt.show()