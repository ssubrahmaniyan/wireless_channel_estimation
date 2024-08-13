# importing the libraries used
import numpy as np, json, math, random, librosa, matplotlib.pyplot as plt, soundfile as sf, sounddevice as sd
from scipy.linalg import lstsq

# initialization of the base variables

Fs = 10000  # Sampling frequency
Fc = 1000000000   # Center frequency
Fd = 100    # Doppler frequency
Tc = 1 / (10*Fd) # Coherence time is 1/4th of the inverse of the doppler frequency

N = 100000  # Number of data points in the signal
M = 3  # Number of different channel values
dom = [-1, 1]
packet_size = math.floor(Fs*Tc)  # Size of each packet
noise_level = 0.001  # noise level(this is 1/SNR in number, the decibel value is 10*log10(SNR))
header_packet_size = math.floor((1/8)*packet_size) # The part of the total signal that is agreed upon and known

# Generate the random signal
def generate_signal(N, mean, var):
    p = (mean + 1) / 2  # Probability of +1, fixes the proportion of 1s that should be spawned in the end
    return np.random.choice([-1, 1], size=N, p=[1 - p, p])

# function to play the different forms of the audio when needed, the argument is the channel index
def playaudio(i, sr):
    print("Playing scaled")
    sd.play(scaled_results[i], sr)
    sd.wait()
    print("playing noisy scaled")
    sd.play(noisy_results[i], sr)
    sd.wait()
    print("playing cleaned")
    sd.play(cleaned_results[i], sr)
    sd.wait()

# function to digitize the values: by using a function similar to logistic regression
def digitize(signal):
    output_signal = np.zeros_like(signal)
    for i, data in enumerate(signal):
        if data >= 0:
            output_signal[i] = 1
        else:
            output_signal[i] = -1
    return output_signal

#use these lines when the audio file is being used for analysis or testing      
#audio_path = 'harvard.wav'
#signal, sr = librosa.load(audio_path)

# use this when a random digital signal is to be created with mean and variance specified
signal = generate_signal(N, 0, 1)

# Load channel values
with open("channels.txt", "r") as f:
    channel_values = [json.loads(line) for line in f]

#channel_values = channel_values[200:]

if len(channel_values) < len(signal):
    print("Warning: Not enough channel values. Repeating the available values.")
    channel_values = (channel_values * (len(signal) // len(channel_values) + 1))[:len(signal)]

# Convert to np array
channel_values = np.array(channel_values)
#channel_values_repeated = np.repeat(channel_values, packet_size, axis=0)[:len(signal)]
channel_values_repeated = channel_values
scaled = signal * channel_values_repeated[:, 0]
#average_signal_power = np.mean(scaled**2) #computing the noise using the SNR of the whole signal
#that assumptions means that even when nothing is being transmit, noise will be
#is that a reasonable model, remains to be seen
#var_noise = average_signal_power * noise_level

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

        #use these lines when the noise is computed with respect to each packet, and remove the previous noise creation lines
        average_signal_power = np.mean(scaled**2)
        var_noise = average_signal_power * noise_level
        
        noise = np.random.normal(0, math.sqrt(var_noise), current_signal.shape)
        noisy = scaled + noise
        
        #Performing LS separation using lstsq
        #the lstsq tries to fit Ax = y
        #because y = hx + n, y/A = x + n/A
        #this is a very good approximation when the snr is relatively row
        #a good value of snr seems to be about 30 dB, which is extremely good
        noisy_dig = digitize(noisy)
        coeff, residuals, rank, s = lstsq(current_signal[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size].reshape(-1, 1), noisy_dig[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size])
        coeff = digitize(coeff)
        cleaned_temp = noisy_dig/coeff
        cleaned = digitize(cleaned_temp) #digitizing the final retrieved output using lstsq
        
        for channel in list(current_channel_values[:, i]):
            #print(channel)
            header_error[i] += (coeff - channel)**2
        header_wrong_count[i] += np.sum(current_signal[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size] != cleaned[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size])
        num[i] += len(current_signal[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size])
        header_bit_count[i] += len(current_channel_values[:, i])

        cleaned_results[i][packet_start:packet_end] = cleaned
        scaled_results[i][packet_start:packet_end] = scaled        
        noisy_results[i][packet_start:packet_end] = noisy


#playaudio(0, sr) #use only when the sample analysed is an audio signal, ensure the functions to load the audio is already run

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
#Compute MSE and correlation with respect to the original signal
# Compute the fraction of corrupted bits with respect to the original signal
for i in range(M):
    total_bits = len(signal)
    corrupted_bits = np.sum(signal != cleaned_results[i])
    fraction_corrupted = corrupted_bits / total_bits
    frac = header_wrong_count[i] / num[i]
    print(f'Channel {i+1}: Fraction of Corrupted Bits(full signal) = {fraction_corrupted:.5f}, MSE(channel values) = {(header_error[i]/header_bit_count[i]):.5f}, Fraction of corrupted bits(Pilot) = {frac:.9f}')
