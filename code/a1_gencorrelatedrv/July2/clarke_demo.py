import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft, fftfreq

Fc = 1  # Carrier frequency in Hz
Fd = 0.00001  # Doppler spread frequency in Hz
t = np.linspace(0, 1, 10000)  # Time vector
fs = len(t)  # Sampling frequency

#two iid random variables
rv1 = np.random.normal(0, 1, len(t))
rv2 = np.random.normal(0, 1, len(t))

#the combined signal
combined_signal = rv1 * np.sin(2 * np.pi * Fc * t) + rv2 * np.cos(2 * np.pi * Fc * t)

#filter using the given PSD
def doppler_filter_psd(signal, Fc, Fd, fs):
    N = len(signal)
    freqs = fftfreq(N, d=1/fs)
    H = np.zeros_like(freqs, dtype=complex)
    
    for i, f in enumerate(freqs):
        if np.abs(f - Fc) < Fd:
            H[i] = 1.5 / (np.pi * Fd * np.sqrt(1 - ((f - Fc) / Fd)**2))
        else:
            H[i] = 0
    
    H = np.nan_to_num(H)
    
    signal_fft = fft(signal)
    filtered_signal_fft = signal_fft * H
    filtered_signal = ifft(filtered_signal_fft)
    
    return filtered_signal

filtered_signal = doppler_filter_psd(combined_signal, Fc, Fd, fs)


combined_sqrt_signal = np.sqrt(rv1**2 + rv2**2)

def autocorrelation(signal):
    result = np.correlate(signal, signal, mode='full')
    return result

# Plotting the filtered signal and its properties
plt.figure(figsize=(12, 16))

# Plot the real part of the filtered signal
plt.subplot(4, 2, 1)
plt.plot(t, np.real(filtered_signal), label='Real Part')
plt.title('Real Part of Filtered Signal')
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.legend()

# Plot the imaginary part of the filtered signal
plt.subplot(4, 2, 2)
plt.plot(t, np.imag(filtered_signal), label='Imaginary Part')
plt.title('Imaginary Part of Filtered Signal')
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.legend()

# Plot the PDF of rv1
plt.subplot(4, 2, 3)
plt.hist(rv1, bins=300, density=True, alpha=0.6, color='g', label='rv1')
plt.title('PDF of Random Variable 1')
plt.xlabel('Value')
plt.ylabel('Probability Density')
plt.legend()

# Plot the PDF of rv2
plt.subplot(4, 2, 4)
plt.hist(rv2, bins=300, density=True, alpha=0.6, color='b', label='rv2')
plt.title('PDF of Random Variable 2')
plt.xlabel('Value')
plt.ylabel('Probability Density')
plt.legend()

# Plot the PDF of combined signal
plt.subplot(4, 2, 5)
plt.hist(combined_signal, bins=300, density=True, alpha=0.6, color='purple', label='Combined Signal')
plt.title('PDF of Combined Signal')
plt.xlabel('Value')
plt.ylabel('Probability Density')
plt.legend()

# Plot the PDF of the combined sqrt signal
plt.subplot(4, 2, 6)
plt.hist(combined_sqrt_signal, bins=300, density=True, alpha=0.6, color='orange', label='Combined Sqrt Signal')
plt.title('PDF of Combined Sqrt Signal')
plt.xlabel('Value')
plt.ylabel('Probability Density')
plt.legend()

# Plot the autocorrelation of the combined signal
plt.subplot(4, 2, 7)
autocorr_combined = autocorrelation(combined_signal)
lags = np.arange(-len(combined_signal) + 1, len(combined_signal))
plt.plot(lags, autocorr_combined, label='Autocorrelation of Combined Signal')
plt.title('Autocorrelation of Combined Signal')
plt.xlabel('Lag')
plt.ylabel('Autocorrelation')
plt.legend()

# Plot the autocorrelation of the filtered signal
plt.subplot(4, 2, 8)
autocorr_filtered = autocorrelation(filtered_signal)
lags = np.arange(-len(filtered_signal) + 1, len(filtered_signal))
plt.plot(lags, autocorr_filtered, label='Autocorrelation of Filtered Signal')
plt.title('Autocorrelation of Filtered Signal')
plt.xlabel('Lag')
plt.ylabel('Autocorrelation')
plt.legend()

plt.tight_layout()
plt.show()

