import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from sympy import Matrix
from statsmodels.tsa.arima_process import ArmaProcess
from scipy.signal import correlate, convolve
from scipy.fft import fft, ifft, fftfreq

N = 100  # Total number of sample points per random variable (component of the random vector)
K = 3     # Number of random variables per random process (number of components)
fd = 0.1  # Maximum Doppler frequency for the Jakes spectrum

def generate_positive_definite_matrix(k):
    A = np.random.rand(k, k)
    A = np.dot(A, A.T)  # Make it symmetric and positive definite
    return A

def generate_arma_process(ar_params, ma_params, n):
    ar = np.array(ar_params)
    ma = np.concatenate(([1], ma_params))  # Adding 1 to the MA parameters for the ARMA model
    arma_process = ArmaProcess(ar, ma)
    return arma_process.generate_sample(n)

# Initialize ARMA processes for each component of Y
phi = 0.6
theta = 0.4
arma_processes = [ArmaProcess([1, -phi], [1, theta]) for _ in range(K)]

# Random mean array for A
A = np.random.randn(K)

# Random positive definite covariance matrix Cxx
Cxx = generate_positive_definite_matrix(K)

cov_mat = Matrix(Cxx)
P, D = cov_mat.diagonalize()
P = np.array(P).astype('float')
D = np.array(D).astype('float')
Drt = np.sqrt(D)

B = P @ Drt

# Initialize arrays for streaming data
Y = np.zeros((K, N))
X = np.zeros((K, N))

# Define H(f) for the Jakes spectrum
def H(f, fd):
    l = np.zeros((len(f),))
    for i in range(len(f)):
        if np.abs(f[i]) < fd:
            l[i] = 1 / np.sqrt(1 - (f[i]/fd)**2)
        else:
            l[i] = 0
    return l

# Generate the frequency domain representation of the filter
freqs = fftfreq(N, d=1/N)
H_f = H(freqs, fd)

# Compute the time domain filter h[n] using inverse FFT
h = ifft(H_f).real

# Set up the figure and axes
fig, axs = plt.subplots(2, 1, figsize=(12, 12), gridspec_kw={'height_ratios': [2, 1]})

# Lines for the time series
lines = [axs[0].plot([], [], marker='o', markersize=2, label=f'Component {k+1}')[0] for k in range(K)]
axs[0].set_xlim(0, N)
axs[0].set_ylim(-5, 5)
axs[0].set_title('Time Series')
axs[0].legend()

# Lines for the cross-correlation plots
cross_corr_lines = [axs[1].plot([], [], label=f'Cross-correlation between Component {i+1} and {j+1}')[0] 
                    for i in range(K) for j in range(i+1, K)]
axs[1].set_xlim(-N + 1, N - 1)
axs[1].set_ylim(-100, 100)
axs[1].set_title('Cross-correlations')
axs[1].legend()

def init():
    for line in lines:
        line.set_data([], [])
    for line in cross_corr_lines:
        line.set_data([], [])
    return lines + cross_corr_lines

def update(frame):
    # Generate new sample for each ARMA process
    new_sample = np.array([arma_process.generate_sample(nsample=1) for arma_process in arma_processes]).flatten()
    Y[:, frame] = new_sample
    
    # Update X with new sample
    X[:, frame] = A + B @ Y[:, frame]
    
    # Convolve each component with the Jakes spectrum filter
    X_filtered = np.zeros_like(X)
    for k in range(K):
        X_filtered[k, :] = convolve(X[k, :], h, mode='same')
    
    for k, line in enumerate(lines):
        line.set_data(np.arange(frame + 1), X_filtered[k, :frame + 1])
    
    # Update cross-correlations
    cross_corr_idx = 0
    for i in range(K):
        for j in range(i+1, K):
            if frame > 0:
                lags, correlation = compute_cross_correlation(X_filtered[i, :frame + 1], X_filtered[j, :frame + 1])
                cross_corr_lines[cross_corr_idx].set_data(lags, correlation)
            cross_corr_idx += 1
    
    return lines + cross_corr_lines

def compute_cross_correlation(x, y):
    correlation = correlate(x, y, mode='full')
    lags = np.arange(-len(x) + 1, len(x))
    return lags, correlation

ani = FuncAnimation(fig, update, frames=range(N), init_func=init, blit=True)
plt.tight_layout()
plt.show()
