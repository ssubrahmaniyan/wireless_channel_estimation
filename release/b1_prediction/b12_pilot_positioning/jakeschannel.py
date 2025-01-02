import numpy as np
from sympy import Matrix
import json
import argparse

def generate_positive_definite_matrix(k):
    A = np.random.rand(k, k)
    A = np.dot(A, A.T)  # Make it symmetric and positive definite
    return A

def jakes_sos(P, K, Fs, Fd, N_paths, typ):
    """
    Modified Jakes simulator implementation that preserves the original function signature
    while implementing equations (7) and (8) from the paper.
    
    Parameters:
    P : int
        Number of time points
    K : int
        Number of realizations
    Fs : float
        Sampling frequency
    Fd : float
        Doppler frequency
    N_paths : int
        Number of paths (will be used to determine M where N = 4M + 2)
    typ : str
        'comp' for complex output, anything else for real output
        
    Returns:
    jakes_rvs : complex array of shape (K, P)
        Array of K realizations of the fading process
    """
    # Generate time vector
    t = np.linspace(0, P/Fs, P)
    wd = 2 * np.pi * Fd
    
    # Calculate M from N_paths (approximating N ≈ N_paths)
    M = (N_paths - 2) // 4
    N = 4 * M + 2  # Actual N used in equations
    
    # Initialize output array
    jakes_rvs = np.zeros((K, P), dtype=complex)
    
    # Generate realizations
    for k in range(K):
        # Generate beta_n according to equation (7c)
        beta = np.zeros(M + 1)
        beta[0] = np.pi/4
        beta[1:] = np.pi * np.arange(1, M + 1) / M
        
        # Generate a_n according to equation (7a)
        a = np.zeros(M + 1)
        a[0] = np.sqrt(2) * np.cos(beta[0])
        a[1:] = 2 * np.cos(beta[1:])
        
        # Generate b_n according to equation (7b)
        b = np.zeros(M + 1)
        b[0] = np.sqrt(2) * np.sin(beta[0])
        b[1:] = 2 * np.sin(beta[1:])
        
        # Generate w_n according to equation (7d)
        w = np.zeros(M + 1)
        w[0] = wd
        w[1:] = wd * np.cos(2 * np.pi * np.arange(1, M + 1) / N)
        
        # Generate random phases for each realization
        phi = np.random.uniform(-np.pi, np.pi, M + 1)
        
        # Calculate u_c(t) according to equation (8b)
        u_c = (2/np.sqrt(N)) * np.sum([a[n] * np.cos(w[n]*t + phi[n]) for n in range(M + 1)], axis=0)
        
        if typ == 'comp':
            # Calculate u_s(t) according to equation (8c)
            u_s = (2/np.sqrt(N)) * np.sum([b[n] * np.cos(w[n]*t + phi[n]) for n in range(M + 1)], axis=0)
            jakes_rvs[k] = u_c + 1j * u_s
        else:
            jakes_rvs[k] = u_c + 1j * 0
            
    return jakes_rvs

def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(description="Generate correlated Jakes random variables and save them to a file.")
    
    # Adding arguments

    parser.add_argument("-P", type=int, default = 100000, help="Number of data lines to generate")
    parser.add_argument("-Fd", type=int, default = 100, help="Doppler spread")
    parser.add_argument("-N", type=int, default=100, help="Number of sample points (rolling window size)")
    parser.add_argument("-K", type=int, default=1, help="Number of random variables per random process (number of components)")
    parser.add_argument("-T", type=str, default = 'real', help = "The data type of the channel values are real or complex")
    parser.add_argument("-Fs", type=int, default=100000, help="Sampling frequency of the channel, measures the data rate effectively")
    # Parse arguments
    args = parser.parse_args()
    
    # Assign arguments to variables
    P = args.P
    Fd = args.Fd
    N = args.N
    K = args.K
    typ = args.T
    Fs = args.Fs  # Sampling frequency
    Fc = 1000000000  # Center frequency
    F = 2 * Fd  # Signal frequency band

    # Generate Jakes random variables using SoS method
    jakes_rvs = jakes_sos(P, K, Fs, Fd, N, typ) # Use 'real' for purely real channels and 'comp' for complex channels

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
    #correlated_jakes_rvs = B @ jakes_rvs
    correlated_jakes_rvs = jakes_rvs
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

    #print(f"Generated {P} data lines and saved them to channels.txt")

if __name__ == "__main__":
    main()
