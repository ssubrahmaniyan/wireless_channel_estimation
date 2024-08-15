import numpy as np, json, math, matplotlib.pyplot as plt
from scipy.linalg import lstsq
from sympy import Matrix
import json

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

def gen_correlated(dopp):    
    # Assign arguments to variables
    P = 100000
    Fd = dopp
    N = 1000
    K = 1
    typ = 'real'
    Fs = 100000  # Sampling frequency
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

def transmit(dopp, noise_level):
    gen_correlated(dopp)
    Fs = 10000  # Sampling frequency
    Fc = 1000000000   # Center frequency
    Fd = dopp    # Doppler frequency
    Tc = 1 / (5*Fd) # Coherence time is 1/4th of the inverse of the doppler frequency

    N = 100000  # Number of data points in the signal
    M = 1  # Number of different channel values
    dom = [-1, 1]
    packet_size = math.floor(Fs*Tc)  # Size of each packet
    noise_level = noise_level  # noise level(this is 1/SNR in number, the decibel value is 10*log10(SNR))
    header_packet_size = math.floor((1/6)*packet_size) # The part of the total signal that is agreed upon and known

    # Generate the random signal
    def generate_signal(N, mean, var):
        p = (mean + 1) / 2  # Probability of +1, fixes the proportion of 1s that should be spawned in the end
        return np.random.choice([-1, 1], size=N, p=[1 - p, p])

    # function to digitize the values: by using a function similar to logistic regression
    def digitize(signal):
        output_signal = np.zeros_like(signal)
        for i, data in enumerate(signal):
            if data >= 0:
                output_signal[i] = 1
            else:
                output_signal[i] = -1
        return output_signal

    # use this when a random digital signal is to be created with mean and variance specified
    signal = generate_signal(N, 0, 1)

    with open("channels.txt", "r") as f:
        # Deserialize each line into complex numbers
        channel_values = [
            np.array([re for re, im in json.loads(line)])
            for line in f
        ]

    if len(channel_values) < len(signal):
        print("Warning: Not enough channel values. Repeating the available values.")
        channel_values = (channel_values * (len(signal) // len(channel_values) + 1))[:len(signal)]

    # Convert to np array
    channel_values = np.array(channel_values)
    channel_values_repeated = channel_values

    # Initialize cleaned results for each channel
    cleaned_results = [np.zeros(len(signal)) for _ in range(M)]
    scaled_results = [np.zeros(len(signal), dtype=np.float32) for _ in range(M)]
    noisy_results = [np.zeros(len(signal), dtype=np.float32) for _ in range(M)]

    header_error = np.zeros(M)
    header_bit_count = np.zeros(M)
    header_wrong_count = np.zeros(M)
    num = np.zeros(M)

    # Processing packet-wise
    for packet_start in range(0, len(signal), packet_size):
        packet_end = min(packet_start + packet_size, len(signal))
        
        current_signal = signal[packet_start:packet_end]
        current_channel_values = channel_values_repeated[packet_start:packet_end]
        for i in range(M):    
            scaled = current_signal * current_channel_values[:, i]

            # Noise computation with respect to each packet
            average_signal_power = np.mean(scaled**2)
            var_noise = average_signal_power * noise_level
            
            noise = np.random.normal(0, math.sqrt(var_noise), current_signal.shape)
            noisy = scaled + noise
            
            # Performing LS separation using lstsq
            noisy_dig = digitize(noisy)
            coeff, residuals, rank, s = lstsq(current_signal[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size].reshape(-1, 1), noisy_dig[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size])
            coeff = digitize(coeff)
            cleaned_temp = noisy_dig/coeff
            cleaned = digitize(cleaned_temp) # digitizing the final retrieved output using lstsq
            
            for channel in list(current_channel_values[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size, i]):
                header_error[i] += (coeff - channel)**2
                header_wrong_count[i] += np.sum(current_signal[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size] != cleaned[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size])
                num[i] += len(current_signal[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size])
                header_bit_count[i] += len(current_channel_values[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size, i])


            cleaned_results[i][packet_start:packet_end] = cleaned
            scaled_results[i][packet_start:packet_end] = scaled        
            noisy_results[i][packet_start:packet_end] = noisy

    # Compute MSE and BER for headers alone and data bits alone
    header_mse = np.zeros(M)
    data_ber = np.zeros(M)

    # Compute metrics
    for i in range(M):
        # Header MSE (already computed, but this is more explicit)
        header_mse[i] = header_error[i] / header_bit_count[i]
        
        # Data bits (excluding header)
        total_data_bits = 0
        corrupted_data_bits = 0
        
        for packet_start in range(0, len(signal), packet_size):
            data_start = packet_start + int(packet_size / 2 + header_packet_size / 2)
            data_end = min(packet_start + packet_size, len(signal))
            
            data_bits = signal[data_start:data_end]
            cleaned_data_bits = cleaned_results[i][data_start:data_end]
            
            # Accumulate total bits and corrupted bits
            total_data_bits += len(data_bits)
            corrupted_data_bits += np.sum(data_bits != cleaned_data_bits)
        
        # Compute BER for data bits alone
        data_ber[i] = corrupted_data_bits / total_data_bits

    return header_mse[0], data_ber[0]

doppler_levels = [5, 10, 20, 50, 75, 100, 150, 200]
h_mses = []
d_bers = []
for doppler_level in doppler_levels:
    h_mse, d_ber = transmit(doppler_level, 0.1)
    h_mses.append(h_mse)
    d_bers.append(d_ber)

plt.plot(doppler_levels, h_mses, label = 'Header MSEs', marker = 'x')
plt.plot(doppler_levels, d_bers, label = 'Data BERs', marker = 'o')
plt.legend()
plt.xlabel("Doppler frequency in Hz")
plt.ylabel('Error value')
plt.title("Plot of MSE, BER for fixed SNR of 0.1 as a function of dooppler frequencies in Hz")
plt.show()
