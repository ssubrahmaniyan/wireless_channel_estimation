import numpy as np
import matplotlib.pyplot as plt
from sympy import Matrix, sqrt
from statsmodels.tsa.arima_process import ArmaProcess
from scipy.signal import correlate

N = 1000  # Number of sample points per random variable (component of the random vector)
K = 3     # Number of random variables per random process (number of components)

def generate_positive_definite_matrix(k):
    # Generate a random lower triangular matrix
    L = np.tril(np.random.rand(k, k))

    # Construct the positive definite matrix by multiplying L with its transpose
    A = np.dot(L, L.T)

    return A

# Function to generate ARMA process
def generate_arma_process(ar_params, ma_params, n):
    ar = np.array(ar_params)
    ma = np.concatenate(([1], ma_params))  # Adding 1 to the MA parameters for the ARMA model
    arma_process = ArmaProcess(ar, ma)
    return arma_process.generate_sample(n)

# Generate ARMA processes for each component of Y
Y = np.zeros((K, N))
for k in range(K):
    phi = 0.6
    theta = 0.4
    arma_sample = generate_arma_process([phi], [theta], N)
    Y[k] = arma_sample

# Random mean array for A
A = np.random.randn(K)

# Random positive definite covariance matrix Cxx
#Cxx = generate_positive_definite_matrix(K)
Cxx = np.mat([[1, 0.2, 0.5],[0.2, 1, 0.1], [0.5, 0.1, 1]])

cov_mat = Matrix(Cxx)
P, D = cov_mat.diagonalize()
P = np.mat(P).astype('float')
D = np.mat(D).astype('float')
Drt = np.sqrt(D)

B = P @ Drt
A_matrix = np.tile(A[:, np.newaxis], (1, N))
X = A_matrix + B @ Y

# Plotting
plt.figure(figsize=(12, 6))
for k in range(K):
    plt.plot(np.array(X[k]).flatten(), marker='o', markersize=2, label=f'Component {k+1}')

plt.title('Random Vector Components')
plt.xlabel('Time')
plt.ylabel('Value')
plt.legend()
plt.show()

# Calculate and plot cross-correlation between each pair of components
plt.figure(figsize=(12, 6))
lags = np.arange(-N + 1, N)  # All possible lags

for i in range(K):
    for j in range(i + 1, K):
        ccf = correlate(X[i], X[j], mode='full') / np.sqrt(np.var(X[i]) * np.var(X[j]))
        mid = len(ccf) // 2
        plt.plot(lags, ccf[mid - N + 1: mid + N].T, label=f'CCF of Component {i+1} and Component {j+1}')

plt.title('Cross-Correlation Functions between Components')
plt.xlabel('Lag')
plt.ylabel('CCF')
plt.legend()
plt.grid(True)
plt.show()
