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
noise_level = 0.05  # Noise level
header_packet_sizes = [4, 8, 16, 32, 64, 128]  # Different header packet sizes to test

# Generate the random signal
def generate_signal(N, mean, var):
    p = (mean + 1) / 2  # Probability of +1
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

# Initialize list to store MSE for each header packet size
mse_results = [[] for _ in range(M)]  # One list for each channel

# Processing for different header packet sizes
for header_packet_size in header_packet_sizes:
    cleaned_results = [np.zeros(len(signal), dtype=np.float32) for _ in range(M)]
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
            
            # Performing LS separation using lstsq
            coeff, residuals, rank, s = lstsq(current_signal[:header_packet_size].reshape(-1, 1), noisy[:header_packet_size])
            cleaned = coeff * current_signal
            
            cleaned_results[i][packet_start:packet_end] = cleaned
            scaled_results[i][packet_start:packet_end] = scaled

    # Compute MSE for each channel
    for i in range(M):
        mse = np.mean((scaled_results[i] - cleaned_results[i])**2)
        mse_results[i].append(mse)

# Plotting MSE vs. Header Packet Size for all channels
plt.figure(figsize=(8, 5))
for i in range(M):
    plt.plot(header_packet_sizes, mse_results[i], marker='o', label=f'Channel {i + 1}')

plt.title('MSE vs. Header Packet Size for All Channels')
plt.xlabel('Header Packet Size')
plt.ylabel('Mean Squared Error (MSE)')
plt.xticks(header_packet_sizes)
plt.grid()
plt.legend()
plt.show()

# Print results
for header_packet_size, mse in zip(header_packet_sizes, mse_results):
    for i in range(M):
        print(f'Header Packet Size: {header_packet_size}, Channel {i + 1} MSE: {mse[i]:.4f}')
