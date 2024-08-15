# importing the libraries used
import numpy as np, json, math, matplotlib.pyplot as plt
from scipy.linalg import lstsq

# initialization of the base variables
Fs = 10000  # Sampling frequency
Fc = 1000000000   # Center frequency
Fd = 100    # Doppler frequency
Tc = 1 / (7*Fd) # Coherence time is 1/4th of the inverse of the doppler frequency

N = 100000  # Number of data points in the signal
M = 3  # Number of different channel values
dom = [-1, 1]
packet_size = math.floor(Fs*Tc)  # Size of each packet
noise_level = 0.01  # noise level(this is 1/SNR in number, the decibel value is 10*log10(SNR))
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
        
        for channel in list(current_channel_values[:, i]):
            header_error[i] += (coeff - channel)**2
        header_wrong_count[i] += np.sum(current_signal[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size] != cleaned[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size])
        num[i] += len(current_signal[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size])
        header_bit_count[i] += len(current_channel_values[:, i])

        cleaned_results[i][packet_start:packet_end] = cleaned
        scaled_results[i][packet_start:packet_end] = scaled        
        noisy_results[i][packet_start:packet_end] = noisy

# Plotting the cleaned results against the original signal
plt.figure(figsize=(12, 8))
for i in range(M):
    plt.subplot(M, 1, i + 1)
    plt.plot(signal[10000:11100] - cleaned_results[i][10000:11100], label='Residuals', alpha=0.5)
    plt.title(f'Comparison of Original and Cleaned Signal (Channel {i+1})')
    plt.legend()
    plt.xlabel('Sample Index')
    plt.ylabel('Amplitude')

plt.tight_layout()
plt.show()

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

# Display the results
for i in range(M):
    total_bits = len(signal)
    corrupted_bits = np.sum(signal != cleaned_results[i])
    fraction_corrupted = corrupted_bits / total_bits
    frac = header_wrong_count[i] / num[i]
    
    print(f'Channel {i+1}: Fraction of Corrupted Bits(full signal) = {fraction_corrupted:.5f}, MSE(channel values) = {(header_error[i]/header_bit_count[i]):.5f}, Fraction of corrupted bits(Pilot) = {frac:.9f}')
    print(f'Channel {i+1}: Header MSE = {header_mse[i]:.5f}, Data BER = {data_ber[i]:.9f}')
