import numpy as np
from sympy import Matrix
import json
import sys

N = 1000  # Number of sample points to display (rolling window size)
K = 3    # Number of random variables per random process (number of components)
Fs = 10000  # Sampling frequency
Fc = 1000000000  # Center frequency
Fd = 100  # Doppler spread
F = 200  # Signal frequency band
P = int(sys.argv[1])  # Number of data lines to generate

def generate_positive_definite_matrix(k):
    A = np.random.rand(k, k)
    A = np.dot(A, A.T)  # Make it symmetric and positive definite
    return A

def jakes_sos_real(P, K, Fs, Fd, N):
    t = np.linspace(0, P / Fs, P)
    omega_d = 2 * np.pi * Fd    
    jakes_rvs = np.zeros((K, P))
    for k in range(K):
        alpha_m = np.random.uniform(0, 2 * np.pi, N)
        a_m = np.random.uniform(0, 2*np.pi, N)
        cos_terms = np.cos(omega_d * t[:, None] * np.cos(alpha_m) + a_m)
        print(cos_terms.shape)
        real_part = np.sqrt(1 / N) * np.sum(cos_terms, axis=1)
        print(real_part.shape)
        jakes_rvs[k] = real_part
    
    return jakes_rvs

# Generate Jakes random variables using SoS method (real part only)
jakes_rvs = jakes_sos_real(P, K, Fs, Fd, N)

# Generate a random positive definite covariance matrix Cxx
Cxx = generate_positive_definite_matrix(K)

# Perform eigenvalue decomposition
cov_mat = Matrix(Cxx)
P_mat, D = cov_mat.diagonalize()
P_mat = np.array(P_mat).astype('float')
D = np.array(D).astype('float')
Drt = np.sqrt(D)

B = P_mat @ Drt

# Induce correlation between Jakes random variables
correlated_jakes_rvs = B @ jakes_rvs

# Open file for writing
with open("channels.txt", "w") as fo:
    for i in range(P):
        # Take the i-th column as the new sample
        new_sample = correlated_jakes_rvs[:, i]
        
        # Write the correlated Jakes random variables to the file
        json.dump(new_sample.tolist(), fo)
        fo.write('\n')  # Add a newline after each data point

print(f"Generated {P} data lines and saved them to channels.txt")
