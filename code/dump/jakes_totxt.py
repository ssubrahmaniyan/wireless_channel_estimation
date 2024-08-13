import numpy as np
from sympy import Matrix
from statsmodels.tsa.arima_process import ArmaProcess
import json
import matplotlib.pyplot as plt
import sys

N = 100  # Number of sample points to display (rolling window size)
K = 3    # Number of random variables per random process (number of components)
Fs = 10000  # Sampling frequency
Fc = 1000000000   # Center frequency
Fd = 100    # Doppler spread
F = 200 # Signal frequency band
P = 0  # Number of data lines to generate

X = np.zeros((K, N))
X_filtered = np.zeros((K, N))

def generate_positive_definite_matrix(k):
    A = np.random.rand(k, k)
    A = np.dot(A, A.T)  # Make it symmetric and positive definite
    return A

def jakes_psd(f, fc, fd):
    # Initialize the PSD array
    psd = np.zeros_like(f)

    # Create a mask for frequencies within the Doppler spread
    mask = np.abs(f - fc) <= fd

    # Calculate PSD for the frequencies within the mask
    psd[mask] = 1 / (np.pi * fd * np.sqrt(1 - ((f[mask] - fc) / fd) ** 2))

    # Also apply to the negative frequency counterpart
    negative_mask = np.abs(f + fc) <= fd
    psd[negative_mask] = 1 / (np.pi * fd * np.sqrt(1 - ((f[negative_mask] + fc) / fd) ** 2))

    # Handle potential division by zero (resulting in NaNs)
    psd[np.isnan(psd)] = 0

    return psd

def plot_filter(h, Fs):
    plt.figure(figsize=(12, 8))

    # Time domain
    plt.subplot(2, 1, 1)
    t = np.arange(len(h)) / Fs
    plt.plot(t, h)
    plt.title('Filter Impulse Response')
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')

    # Frequency domain
    plt.subplot(2, 1, 2)
    f = np.linspace(-Fs/2, Fs/2, len(h))
    H = np.fft.fftshift(np.fft.fft(h))
    plt.plot(f, np.abs(H))
    plt.title('Filter Frequency Response')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude')
    plt.xlim(-Fs/2, Fs/2)

    plt.tight_layout()
    plt.show()

# Initialize ARMA processes for each component of Y
phi = np.random.rand()
theta = np.random.rand()
arma_processes = [ArmaProcess([1, -phi], [1, theta]) for _ in range(K)]
# Random mean array for A
A = np.random.randn(K)

# Random positive definite covariance matrix Cxx
Cxx = generate_positive_definite_matrix(K)

cov_mat = Matrix(Cxx)
P_mat, D = cov_mat.diagonalize()
P_mat = np.array(P_mat).astype('float')
D = np.array(D).astype('float')
Drt = np.sqrt(D)

B = P_mat @ Drt

# Generate filter
f_positive = np.linspace(Fc - F, Fc + F, N // 2)
f_negative = np.linspace(-Fc - F, -Fc + F, N // 2)
f = np.concatenate((f_negative, f_positive))
psd = jakes_psd(f, Fc, Fd)
h = np.fft.ifftshift(np.fft.ifft(np.sqrt(psd)).real)

#plot_filter(h, Fs)
# Initialize arrays for data
Y = np.zeros((K, N))
X = np.zeros((K, N))
X_filtered = np.zeros((K, N))

P = int(sys.argv[1])
# Open file for writing
with open("channels.txt", "w") as fo:
    # Generate P data lines
    for _ in range(P):
        # Generate new sample for each ARMA process
        new_sample = np.array([arma_process.generate_sample(nsample=1) for arma_process in arma_processes]).flatten()
        # Shift the data in Y and X arrays to the left by one and add new sample
        Y[:, :-1] = Y[:, 1:]
        Y[:, -1] = new_sample
        X[:, :-1] = X[:, 1:]
        X[:, -1] = A + B @ new_sample
            
        # Apply filter to the last 100 samples
        for k in range(K):
            X_filtered[k, :-1] = X_filtered[k, 1:]  # Shift existing filtered data
            X_filtered[k, -1] =  np.convolve(X[k, -Fs:], h[-Fs:], mode='valid')[-1]
        
        # Write the filtered data to the file
        json.dump(X_filtered[:, -1].tolist(), fo)
        fo.write('\n')  # Add a newline after each data point

print(f"Generated {P} data lines and saved them to channels.txt")
