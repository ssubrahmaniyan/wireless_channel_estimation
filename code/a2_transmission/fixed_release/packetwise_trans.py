# Importing the necessary libraries for animation
import numpy as np, json, math, matplotlib.pyplot as plt, soundfile as sf
#import librosa, random,  sounddevice as sd
from scipy.linalg import lstsq
from matplotlib.animation import FuncAnimation

# Initialization of the base variables
Fs = 10000  # Sampling frequency
Fc = 1000000000   # Center frequency
Fd = 100    # Doppler frequency
Tc = 1 / (Fd) # Coherence time is 1/4th of the inverse of the doppler frequency

N = 100000  # Number of data points in the signal
M = 3  # Number of different channel values
dom = [-1, 1]
packet_size = math.floor(Fs*Tc)  # Size of each packet
noise_level = 0.1  # noise level(this is 1/SNR in number, the decibel value is 10*log10(SNR))
header_packet_size = math.floor((1/6)*packet_size) # The part of the total signal that is agreed upon and known

# Generate the random signal
def generate_signal(N, mean, var):
    p = (mean + 1) / 2  # Probability of +1, fixes the proportion of 1s that should be spawned in the end
    return np.random.choice([-1, 1], size=N, p=[1 - p, p])

# Function to digitize the values
def digitize(signal):
    output_signal = np.zeros_like(signal)
    for i, data in enumerate(signal):
        if data >= 0:
            output_signal[i] = 1
        else:
            output_signal[i] = -1
    return output_signal

# Function to update the animation
def update_plot(packet_idx):
    plt.clf()
    packet_start = packet_idx * packet_size
    packet_end = min(packet_start + packet_size, len(signal))
    
    current_signal = signal[packet_start:packet_end]
    current_channel_values = channel_values_repeated[packet_start:packet_end]
    
    for i in range(1):    
        scaled = current_signal * current_channel_values[:, i]
        average_signal_power = np.mean(scaled**2)
        var_noise = average_signal_power * noise_level
        noise = np.random.normal(0, math.sqrt(var_noise), current_signal.shape)
        noisy = scaled + noise
        noisy_dig = digitize(noisy)
        coeff, residuals, rank, s = lstsq(current_signal[:header_packet_size].reshape(-1, 1), noisy_dig[:header_packet_size])
        cleaned_temp = noisy_dig / coeff
        cleaned = digitize(cleaned_temp)

        plt.subplot(5, 1, 1)
        plt.plot(current_signal, label='Original Signal')
        plt.plot(current_channel_values[:, i], marker = 'x',  label = 'Channel Values')
        plt.legend()
        plt.title(f'Packet {packet_idx + 1}, Channel {i+1}')
        
        plt.subplot(5, 1, 2)
        plt.plot(scaled, label='Scaled Signal')
        plt.legend()
        
        plt.subplot(5, 1, 3)
        plt.plot(noisy, label='Noisy Signal')
        plt.legend()
        
        plt.subplot(5, 1, 4)
        plt.plot(noisy_dig, label='Noisy Digitized Signal')
        plt.plot(cleaned, label='Cleaned Signal', alpha=0.7)
        plt.legend()
        
        plt.subplot(5, 1, 5)
        plt.plot(current_signal - cleaned, label='Difference (Original - Cleaned)')
        plt.legend()
        

    plt.tight_layout()

# Function to handle key press events
def on_key(event):
    global packet_idx
    if event.key == 'enter':
        packet_idx += 1
        if packet_idx >= len(signal) // packet_size:
            packet_idx = 0
        update_plot(packet_idx)
        plt.draw()

# use this when a random digital signal is to be created with mean and variance specified
signal = generate_signal(N, 0, 1)

with open("channels.txt", "r") as f:
    # Deserialize each line into complex numbers
    channel_values = [
        np.array([complex(re, im) for re, im in json.loads(line)])
        for line in f
    ]

channel_values = np.array([[val.real for val in channel] for channel in channel_values])

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

# Variables for controlling the animation
packet_idx = 0

# Create the plot and connect the key press event handler
fig = plt.figure(figsize=(15, 12))
fig.canvas.mpl_connect('key_press_event', on_key)
update_plot(packet_idx)
plt.show()

# Compute the fraction of corrupted bits with respect to the original signal
for i in range(M):
    total_bits = len(signal)
    corrupted_bits = np.sum(signal != cleaned_results[i])
    fraction_corrupted = corrupted_bits / total_bits
    print(f'Channel {i+1}: Fraction of Corrupted Bits = {fraction_corrupted:.5f}')