import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from sympy import Matrix
from scipy.signal import correlate, convolve
from scipy.stats import gaussian_kde
from scipy.fft import fft, fftfreq
from scipy.special import jv

N = 100  # Number of sample points to display (rolling window size)
K = 3   # Number of random variables per random process (number of components)

def generate_positive_definite_matrix(k):
    A = np.random.rand(k, k)
    A = np.dot(A, A.T)  # Make it symmetric and positive definite
    return A

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

# Set up the figure and axes for the animation
fig, axs = plt.subplots(4, 1, figsize=(12, 24), gridspec_kw={'height_ratios': [2, 1, 1, 1]})

# Lines for the time series
lines_100 = [axs[0].plot([], [], marker='o', markersize=2, label=f'Component {k+1}')[0] for k in range(K)]
axs[0].set_xlim(0, 100)
axs[0].set_ylim(-2, 2)
axs[0].set_title('Time Series (Last 100 points)')
axs[0].legend()

# Lines for the cross-correlation plots
cross_corr_lines = [axs[1].plot([], [], label=f'Cross-correlation between Component {i+1} and {j+1}')[0] 
                    for i in range(K) for j in range(i+1, K)]
axs[1].set_xlim(-N + 1, N - 1)
axs[1].set_ylim(-1, 1)
axs[1].set_title('Cross-correlations')
axs[1].legend()

# Lines for the PDF plots
pdf_lines = [axs[2].plot([], [], label=f'PDF of Component {k+1}')[0] for k in range(K)]
axs[2].set_xlim(-5, 5)
axs[2].set_ylim(0, 1)
axs[2].set_title('PDFs')
axs[2].legend()

# Lines for the Fourier transform plots
fourier_lines = [axs[3].plot([], [], label=f'Fourier Transform of Component {k+1}')[0] for k in range(K)]
axs[3].set_xlim(0, N//2)
axs[3].set_ylim(0, 50)
axs[3].set_title('Rolling Fourier Transform')
axs[3].legend()

def init():
    for line in lines_100:
        line.set_data([], [])
    for line in cross_corr_lines:
        line.set_data([], [])
    for line in pdf_lines:
        line.set_data([], [])
    for line in fourier_lines:
        line.set_data([], [])
    return lines_100 + cross_corr_lines + pdf_lines + fourier_lines

def update(frame):
    # Generate new Rayleigh random sample
    new_sample = np.sqrt(np.random.randn(K, 1)**2 + np.random.randn(K, 1)**2).mean(axis=1)
    
    # Shift the data in Y and X arrays to the left by one and add new sample
    Y[:, :-1] = Y[:, 1:]
    Y[:, -1] = new_sample
    X[:, :-1] = X[:, 1:]
    X[:, -1] = A + B @ new_sample
    
    for k, line_100 in enumerate(lines_100):
        line_100.set_data(np.arange(100), X[k, :])
    
    cross_corr_idx = 0
    for i in range(K):
        for j in range(i+1, K):
            if frame > 0:
                lags, correlation = compute_cross_correlation(X[i, :], X[j, :])
                cross_corr_lines[cross_corr_idx].set_data(lags, correlation)
            cross_corr_idx += 1

    for k, line_pdf in enumerate(pdf_lines):
        if frame > 0:
            kde = gaussian_kde(X[k, :])
            x_grid = np.linspace(-5, 5, 500)
            line_pdf.set_data(x_grid, kde(x_grid))
    
    for k, line_fourier in enumerate(fourier_lines):
        if frame > 0:
            fft_values = np.abs(fft(X[k, :]))[:N//2]
            line_fourier.set_data(np.arange(N//2), fft_values)
    
    return lines_100 + cross_corr_lines + pdf_lines + fourier_lines

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

