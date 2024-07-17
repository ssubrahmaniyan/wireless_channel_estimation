import numpy as np
from scipy.linalg import lstsq
import json
import math
import random
import matplotlib.pyplot as plt

# Configuration
N = 100000  # Number of data points in the signal
M = 3  # Number of different channel values
dom = [-1, 1]
packet_size = 256  # Size of each packet
noise_level = 0.01  # Noise level
header_packet_size = 16 # The part of the total signal that is agreed upon and known

# Generate the random signal
def generate_signal(N, mean, var):
    # Calculate the proportion of +1 and -1 needed
    p = (mean + 1) / 2  # Probability of +1
    # Generate signal
    return np.random.choice([-1, 1], size=N, p=[1 - p, p])

signal = generate_signal(N, 0, 1)

# Load channel values
with open("channels.txt", "r") as f:
    channel_values = [json.loads(line) for line in f]

if len(channel_values) < len(signal):
    print("Warning: Not enough channel values. Repeating the available values.")
    channel_values = (channel_values * (len(signal) // len(channel_values) + 1))[:len(signal)]

# Convert to np array
channel_values = np.array(channel_values)
channel_values_repeated = np.repeat(channel_values, packet_size, axis=0)[:len(signal)]

# Initialize cleaned results for each channel
cleaned_results = [np.zeros(len(signal)) for _ in range(M)]
scaled_results = [np.zeros(len(signal), dtype=np.float32) for _ in range(M)]

# Processing packet-wise
for packet_start in range(0, len(signal), packet_size):
    packet_end = min(packet_start + packet_size, len(signal))
    
    current_signal = signal[packet_start:packet_end]
    current_channel_values = channel_values_repeated[packet_start:packet_end]
    
    for i in range(M):    
        scaled = current_signal * current_channel_values[:, i]
        
        average_signal_power = np.mean(scaled**2)
        var_noise = average_signal_power * noise_level
        
        noise = np.random.normal(0, math.sqrt(var_noise), current_signal.shape)
        noisy = scaled + noise
        
        #Performing LS separation using lstsq
        coeff, residuals, rank, s = lstsq(current_signal[:header_packet_size].reshape(-1, 1), noisy[:header_packet_size])
        cleaned = noisy/coeff
        
        cleaned_results[i][packet_start:packet_end] = cleaned
        scaled_results[i][packet_start:packet_end] = scaled

# Plotting the cleaned results against the original signal
plt.figure(figsize=(12, 8))
for i in range(M):
    plt.subplot(M, 1, i + 1)
    plt.plot(scaled_results[i], label='Original Signal', alpha=0.5)
    plt.plot(cleaned_results[i], label=f'Cleaned Signal (Channel {i+1})', alpha=0.75)
    plt.title(f'Comparison of Original and Cleaned Signal (Channel {i+1})')
    plt.legend()
    plt.xlabel('Sample Index')
    plt.ylabel('Amplitude')

plt.tight_layout()
plt.show()

# Compute MSE and correlation
for i in range(M):
    mse = np.mean((scaled_results[i] - cleaned_results[i])**2)
    correlation = np.corrcoef(signal, cleaned_results[i])[0, 1]
    
    print(f'Channel {i+1}: MSE = {mse:.5f}, Correlation = {correlation:.4f}')
