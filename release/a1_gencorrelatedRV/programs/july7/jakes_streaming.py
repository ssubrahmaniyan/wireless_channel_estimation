import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from sympy import Matrix
from statsmodels.tsa.arima_process import ArmaProcess
from scipy.signal import correlate
from scipy.stats import gaussian_kde
import json

N = 100  # Number of sample points to display (rolling window size)
K = 3    # Number of random variables per random process (number of components)
Fs = 100  # Sampling frequency
Fc = 10   # Center frequency
Fd = 1    # Doppler spread
global X, X_filtered
X = np.zeros((K, N))
X_filtered = np.zeros((K, N))

#fo = open(r"channels", "w")

def generate_positive_definite_matrix(k):
    A = np.random.rand(k, k)
    A = np.dot(A, A.T)  # Make it symmetric and positive definite
    return A

def generate_arma_process(ar_params, ma_params, n):
    ar = np.array(ar_params)
    ma = np.concatenate(([1], ma_params))  # Adding 1 to the MA parameters for the ARMA model
    arma_process = ArmaProcess(ar, ma)
    return arma_process.generate_sample(n)

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

# Generate filter
def jakes_psd(f, fc, fd):
    mask = np.abs(f - fc) <= fd
    psd = np.zeros_like(f)
    psd[mask] = 1 / (np.pi * fd * np.sqrt(1 - ((f[mask] - fc) / fd) ** 2))
    return psd

f = np.linspace(-Fs/2, Fs/2, N)
psd = jakes_psd(f, Fc, Fd)
h = np.fft.ifftshift(np.fft.ifft(np.sqrt(psd)).real)
print(h)

# Plot the filter
plot_filter(h, Fs)

# Initialize arrays for streaming data
Y = np.zeros((K, N))
X = np.zeros((K, N))
X_filtered = np.zeros((K, N))

# Set up the figure and axes
fig, axs = plt.subplots(3, 2, figsize=(15, 20))

# Lines for the time series
lines_100 = [axs[0, 0].plot([], [], marker='o', markersize=2, label=f'Component {k+1}')[0] for k in range(K)]
lines_100_filtered = [axs[0, 1].plot([], [], marker='o', markersize=2, label=f'Component {k+1}')[0] for k in range(K)]
for ax in axs[0]:
    ax.set_xlim(0, 100)
    ax.set_ylim(-5, 5)
    ax.legend()
axs[0, 0].set_title('Time Series (Before Filtering)')
axs[0, 1].set_title('Time Series (After Filtering)')

# Lines for the cross-correlation plots
cross_corr_lines = [axs[1, 0].plot([], [], label=f'Corr {i+1}-{j+1}')[0] 
                    for i in range(K) for j in range(i+1, K)]
cross_corr_lines_filtered = [axs[1, 1].plot([], [], label=f'Corr {i+1}-{j+1}')[0] 
                             for i in range(K) for j in range(i+1, K)]
for ax in axs[1]:
    ax.set_xlim(-N + 1, N - 1)
    ax.set_ylim(-1, 1)
    ax.legend()
axs[1, 0].set_title('Cross-correlations (Before Filtering)')
axs[1, 1].set_title('Cross-correlations (After Filtering)')

# Lines for the PDF plots
pdf_lines = [axs[2, 0].plot([], [], label=f'PDF of Component {k+1}')[0] for k in range(K)]
pdf_lines_filtered = [axs[2, 1].plot([], [], label=f'PDF of Component {k+1}')[0] for k in range(K)]
for ax in axs[2]:
    ax.set_xlim(-5, 5)
    ax.set_ylim(0, 1)
    ax.legend()
axs[2, 0].set_title('PDFs (Before Filtering)')
axs[2, 1].set_title('PDFs (After Filtering)')

def init():
    for line in lines_100 + lines_100_filtered + cross_corr_lines + cross_corr_lines_filtered + pdf_lines + pdf_lines_filtered:
        line.set_data([], [])
    return lines_100 + lines_100_filtered + cross_corr_lines + cross_corr_lines_filtered + pdf_lines + pdf_lines_filtered

def update(frame):
    global X, X_filtered  # Make these global so we can modify them

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
        X_filtered[k, -1] = 10 * np.convolve(X[k, -100:], h[-100:], mode='valid')[-1]
    
#    json.dump(X_filtered[:, -1].tolist(), fo)
            
    for k, (line_100, line_100_filtered) in enumerate(zip(lines_100, lines_100_filtered)):
        line_100.set_data(np.arange(100), X[k])
        line_100_filtered.set_data(np.arange(100), X_filtered[k])

    cross_corr_idx = 0
    for i in range(K):
        for j in range(i+1, K): 
            if frame > 0:
                lags, correlation = compute_cross_correlation(X[i], X[j])
                cross_corr_lines[cross_corr_idx].set_data(lags, correlation)
                lags, correlation = compute_cross_correlation(X_filtered[i], X_filtered[j])
                cross_corr_lines_filtered[cross_corr_idx].set_data(lags, correlation)
            cross_corr_idx += 1

    for k, (line_pdf, line_pdf_filtered) in enumerate(zip(pdf_lines, pdf_lines_filtered)):
        if frame > 0:
            kde = gaussian_kde(X[k])
            x_grid = np.linspace(-5, 5, 500)
            line_pdf.set_data(x_grid, kde(x_grid))
            kde = gaussian_kde(X_filtered[k])
            line_pdf_filtered.set_data(x_grid, kde(x_grid))
    
    return lines_100 + lines_100_filtered + cross_corr_lines + cross_corr_lines_filtered + pdf_lines + pdf_lines_filtered
    
def compute_cross_correlation(x, y):
    x = (x - np.mean(x))/np.std(x)
    y = (y - np.mean(y))/np.std(y)
    correlation = correlate(x, y, mode='full')
    lags = np.arange(-len(x) + 1, len(x))
    correlation = correlation / (len(x))
    return lags, correlation

ani = FuncAnimation(fig, update, frames=np.arange(0, float(100)), init_func=init, blit=True, interval=100)
plt.tight_layout()
plt.show()
#fo.close()
