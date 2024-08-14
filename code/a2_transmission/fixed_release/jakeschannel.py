import numpy as np
from sympy import Matrix
import json
import argparse

def generate_positive_definite_matrix(k):
    A = np.random.rand(k, k)
    A = np.dot(A, A.T)  # Make it symmetric and positive definite
    return A

def jakes_sos(P, K, Fs, Fd, N, typ):
    t = np.linspace(0, P/Fs, P)
    omega_d = 2 * np.pi * Fd
    
    # Initialize jakes_rvs to store real or complex numbers
    jakes_rvs = np.zeros((K, P), dtype=complex)
    
    for k in range(K):
        alpha_m = np.random.uniform(0, 2 * np.pi, N)
        a_m = np.random.uniform(0, 2 * np.pi, N)
        b_m = np.random.uniform(0, 2 * np.pi, N)
        
        cosine_terms = np.cos((omega_d * t[:, None] * np.cos(alpha_m)) + a_m)
        real_part = np.sqrt(1/N) * np.sum(cosine_terms, axis=1)
        
        if typ == 'comp':
            sine_terms = np.sin((omega_d * t[:, None] * np.cos(alpha_m)) + b_m)
            imag_part = np.sqrt(1/N) * np.sum(sine_terms, axis=1)
            jakes_rvs[k] = real_part + 1j * imag_part
        else:
            jakes_rvs[k] = real_part + 1j * 0 # Enforcing the real numbers also to be modelled as complex for easy use
    
    return jakes_rvs

def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(description="Generate correlated Jakes random variables and save them to a file.")
    
    # Adding arguments
    parser.add_argument("-P", type=int, default = 100000, help="Number of data lines to generate")
    parser.add_argument("-Fd", type=int, default = 100, help="Doppler spread")
    parser.add_argument("-N", type=int, default=1000, help="Number of sample points (rolling window size)")
    parser.add_argument("-K", type=int, default=3, help="Number of random variables per random process (number of components)")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Assign arguments to variables
    P = args.P
    Fd = args.Fd
    N = args.N
    K = args.K
    Fs = 100000  # Sampling frequency
    Fc = 1000000000  # Center frequency
    F = 2 * Fd  # Signal frequency band

    # Generate Jakes random variables using SoS method
    jakes_rvs = jakes_sos(P, K, Fs, Fd, N, 'comp') # Use 'real' for purely real channels and 'comp' for complex channels

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

    # Write the data to a file
    with open("channels.txt", "w") as fo:
        for i in range(P):
            # Take the i-th column as the new sample
            new_sample = correlated_jakes_rvs[:, i]
            
            # Convert complex numbers to a format that can be serialized (e.g., tuple of real and imaginary parts)
            new_sample_serialized = [(z.real, z.imag) for z in new_sample]
            
            # Write the serialized complex numbers to the file
            json.dump(new_sample_serialized, fo)
            fo.write('\n')  # Add a newline after each data point

    print(f"Generated {P} data lines and saved them to channels.txt")

if __name__ == "__main__":
    main()
