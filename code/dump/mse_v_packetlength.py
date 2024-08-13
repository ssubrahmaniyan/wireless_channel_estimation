# importing the libraries used
import numpy as np, json, math, random, librosa, matplotlib.pyplot as plt, soundfile as sf, sounddevice as sd, gc
from scipy.linalg import lstsq

# initialization of the base variables
Fs = 10000  # Sampling frequency
Fc = 1000000000   # Center frequency
Fd = 200    # Doppler frequency
Tc = 1 / (4*Fd) # Coherence time is 1/4th of the inverse of the doppler frequency

N = 100000  # Number of data points in the signal
M = 3  # Number of different channel values
dom = [-1, 1]
packet_size_factor = math.floor(Fs*Tc)  
packet_size_ratios = [0.25, 0.5, 1, 2]  
packet_sizes = [math.floor(packet_size_factor * i) for i in packet_size_ratios]
header_packet_size = math.floor((1/8) * packet_size_factor)
noise_level = 0.001

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

def plotter():
    plt.figure(figsize=(8, 5))
    for i in range(M):
        plt.plot(packet_sizes, msearr[:][i], marker='o', label=f'Channel {i + 1}')
    plt.title('MSE vs. Data Packet Size for All Channels With Fixed Header Size/ Packet Size Ratio(1/8)')
    plt.xlabel('Data Packet Size')
    plt.ylabel('Mean Squared Error (MSE)')
    plt.xticks(packet_sizes)
    plt.grid()
    plt.legend()
    plt.show()

# use this when a random digital signal is to be created with mean and variance specified
signal = generate_signal(N, 0, 1)

# Load channel values
with open("channels.txt", "r") as f:
    channel_values = [json.loads(line) for line in f]

if len(channel_values) < len(signal):
    print("Warning: Not enough channel values. Repeating the available values.")
    channel_values = (channel_values * (len(signal) // len(channel_values) + 1))[:len(signal)]

#computing the noise using the SNR of the whole signal
#that assumptions means that even when nothing is being transmit, noise will be
#is that a reasonable model, remains to be seen

msearr = [[] for _ in range(M)]

for i, packet_size in enumerate(packet_sizes):
    # Convert to np array
    channel_values = np.array(channel_values)
    channel_values_repeated = np.repeat(channel_values, packet_size, axis=0)[:len(signal)]

    average_signal_powers = [0 for i in range(M)]

    for i in range(M):
        scaled = signal * channel_values_repeated[:, i]
        average_signal_powers[i] = np.mean(scaled**2)
    
    # Initialize cleaned results for each channel
    cleaned_results = [np.zeros(len(signal)) for _ in range(M)]
    scaled_results = [np.zeros(len(signal), dtype=np.float32) for _ in range(M)]
    noisy_results = [np.zeros(len(signal), dtype=np.float32) for _ in range(M)]

    # Processing packet-wise
    for packet_start in range(0, len(signal), packet_size):
        packet_end = min(packet_start + packet_size, len(signal))
        
        current_signal = signal[packet_start:packet_end]
        current_channel_values = channel_values_repeated[packet_start:packet_end]
        
        for i in range(M):    
            scaled = current_signal * current_channel_values[:, i]

            #use these lines when the noise is computed with respect to each packet, and remove the previous noise creation lines
            #average_signal_power = np.mean(scaled**2)
            #var_noise = average_signal_power * noise_level
            
            noise = np.random.normal(0, math.sqrt(average_signal_powers[i] * noise_level), current_signal.shape)
            noisy = scaled + noise
            
            #Performing LS separation using lstsq
            #the lstsq tries to fit Ax = y
            #because y = hx + n, y/A = x + n/A
            #this is a very good approximation when the snr is relatively row
            #a good value of snr seems to be about 30 dB, which is extremely good
            noisy_dig = digitize(noisy)
            coeff, residuals, rank, s = lstsq(current_signal[:header_packet_size].reshape(-1, 1), noisy_dig[:header_packet_size])
            cleaned_temp = noisy_dig/coeff
            
            cleaned = digitize(cleaned_temp) #digitizing the final retrieved output using lstsq

            cleaned_results[i][packet_start:packet_end] = cleaned
            scaled_results[i][packet_start:packet_end] = scaled        
            noisy_results[i][packet_start:packet_end] = noisy

    #Compute MSE and correlation with respect to the original signal
    for i in range(M):
        mse = np.mean((signal - cleaned_results[i])**2)
        correlation = np.corrcoef(signal, cleaned_results[i])[0, 1]
        print(f'Channel {i+1}: MSE = {mse:.5f}, Correlation = {correlation:.4f}')
        msearr[i].append(mse)
    
    del channel_values_repeated
    del cleaned_results
    del scaled_results
    del noisy_results
    gc.collect()


#print(msearr)
plotter()