import numpy as np, json, math, matplotlib.pyplot as plt
from scipy.linalg import lstsq
from scipy.special import erfc
def generate_signal(N, mean, var):
    p = (mean + 1) / 2  # Probability of +1, fixes the proportion of 1s that should be spawned in the end
    return np.random.choice([-1, 1], size=N, p=[1 - p, p])

with open("channels.txt", "r") as f:
    # Deserialize each line into complex numbers
    channel_vals = [json.loads(line)[0] for line in f]
    channel_vals = np.array([complex(re, im) for re, im in channel_vals])

N = len(channel_vals)
#signal = generate_signal(100000, 0, 1)
ip = np.random.rand(N) > 0.5
signal = 2 * ip - 1
noise_levels = np.arange(-3, 36, 5)

def transmit(dopp, noise_level):
    Fs = 100000  # Sampling frequency
    Fc = 1000000000   # Center frequency
    Fd = dopp    # Doppler frequency
    Tc = 1 / (10*Fd) # Coherence time is 1/4th of the inverse of the doppler frequency
     # Number of data points in the signal
    M = 1  # Number of different channel values
    dom = [-1, 1]
    #packet_size = math.floor(Fs*Tc)  # Size of each packet
    packet_size = 1
    noise_level = noise_level  # noise level(this is 1/SNR in number, the decibel value is 10*log10(SNR))
    #header_packet_size = math.floor((1/4)*packet_size) # The part of the total signal that is agreed upon and known
    header_packet_size = 1
    # Generate the random signal


    # function to digitize the values: by using a function similar to logistic regression
    def digitize(signal):
        output_signal = np.zeros_like(signal)
        for i, data in enumerate(signal):
            if data >= 0:
                output_signal[i] = 1
            else:
                output_signal[i] = 0
        return output_signal

    # use this when a random digital signal is to be created with mean and variance specified

    # Convert to np array
    channel_values_repeated = channel_vals

    # Initialize cleaned results for each channel
    cleaned_results = [np.zeros(len(signal), dtype = np.complex128) for _ in range(M)]
    scaled_results = [np.zeros(len(signal), dtype=np.complex128) for _ in range(M)]
    noisy_results = [np.zeros(len(signal), dtype=np.complex128) for _ in range(M)]

    header_error = np.zeros(M)
    header_bit_count = np.zeros(M)
    header_wrong_count = np.zeros(M)
    num = np.zeros(M)
    
    # noise = (1/np.sqrt(2)) * (np.random.randn(N) + 1j * np.random.randn(N))
    # #h = (1/np.sqrt(2)) * (np.random.randn(N) + 1j * np.random.randn(N))
    # h = np.array(channel_vals)
    # #h = np.zeros_like(channel_vals)
    # #h = [1 for el in h]
    # y = h * signal + 10**(-noise_level / 20) * noise
    # yHat = y / h
    # ipHat = np.real(yHat) > 0
    # nErr = np.sum(ip != ipHat)

    # simBer = nErr
    #Processing packet-wise
    for packet_start in range(0, len(signal), packet_size):
        packet_end = min(packet_start + packet_size, len(signal))
        
        current_signal = signal[packet_start:packet_end]
        current_channel_values = channel_values_repeated[packet_start:packet_end]
        for i in range(M):    
            scaled = current_signal * current_channel_values[:]

            # Noise computation with respect to each packet
            average_signal_power = np.mean((abs(scaled))**2)
            #average_signal_power = np.mean(current_signal**2)
            var_noise = average_signal_power / noise_level

            #noise = np.random.normal(0, np.sqrt(var_noise), current_signal.shape)            
            noise = (10 ** (-noise_level/20))*np.random.randn(1)
            #noisy = scaled + noise
            noisecomp = (10 ** (-noise_level/20))*np.random.randn(1)
            noisy = scaled + np.sqrt(0.5) * complex(noise, noisecomp)
            # Performing LS separation using lstsq
            noisy_dig = noisy
            #noisy_dig = digitize(noisy)
            #coeff, residuals, rank, s = lstsq(current_signal[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size].reshape(-1, 1), noisy_dig[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size])
            coeff = current_channel_values[:]
            #coeff = digitize(coeff)
            # noisy_dig = digitize(noisy_dig)
            cleaned_temp = (noisy_dig/coeff)
            cleaned = digitize(np.real(cleaned_temp)) # digitizing the final retrieved output using lstsq
            
            for channel in list(current_channel_values[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size]):
                header_error[i] += (abs(coeff - channel))**2
                header_wrong_count[i] += np.sum(current_signal[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size] != cleaned[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size])
                num[i] += len(current_signal[(int)(packet_size/2 - header_packet_size / 2):(int)(packet_size/2 - header_packet_size / 2) + header_packet_size])
                header_bit_count[i] += 1


            cleaned_results[i][packet_start:packet_end] = cleaned
            scaled_results[i][packet_start:packet_end] = scaled        
            noisy_results[i][packet_start:packet_end] = noisy
    # plt.plot(noisy_results[0][1000:1050], label = 'noisy')
    # plt.plot(scaled_results[0][1000:1050], label = 'scaled')
    # plt.plot(cleaned_results[0][1000:1050], label = 'cleaned')
    # plt.plot(signal[1000:1050], label = "original")
    # plt.legend()
    # plt.show()
    # Compute MSE and BER for headers alone and data bits alone
    header_mse = np.zeros(M)
    data_ber = np.zeros(M, dtype = np.float128)

    # Compute metrics
    for j in range(M):
        # Header MSE (already computed, but this is more explicit)
        header_mse[j] = header_error[j] / header_bit_count[j]
        
        # Data bits (excluding header)
        total_data_bits = 0
        corrupted_data_bits = 0
        
        for packet_start in range(0, len(signal), packet_size):
            #data_start = packet_start + int(packet_size / 2 + header_packet_size / 2)
            #data_end = min(packet_start + packet_size, len(signal))
            
            data_start1 = packet_start 
            data_end1 = packet_start + (int)(packet_size/2 - header_packet_size / 2) + 1
            #data_start2 = packet_start + (int)(packet_size/2 - header_packet_size / 2) + header_packet_size
            #data_end2 = min(packet_start + packet_size, len(signal))
            #data_bits = np.concatenate((signal[data_start1:data_end1], signal[data_start2:data_end2]), axis = 0)
            #cleaned_data_bits = np.concatenate((cleaned_results[j][data_start1:data_end1], cleaned_results[j][data_start2:data_end2]), axis = 0)
            data_bits = ip[data_start1:data_end1]
            cleaned_data_bits = cleaned_results[j][data_start1:data_end1]
            # Accumulate total bits and corrupted bits
            total_data_bits += len(data_bits)
            corrupted_data_bits += np.sum(data_bits != cleaned_data_bits)
        
        # Compute BER for data bits alone
        # print(signal[0:1000] == cleaned_results[0][0:1000])
        data_ber[j] = np.float128(corrupted_data_bits) / np.float128(total_data_bits)

    return data_ber[0]

    #return simBer

# Desired dB values: 1, 2, 3, ..., 15
#dB_values = np.arange(1, 30, 3)
#noise_levels = dB_values
# Calculate the corresponding noise levels
#noise_levels = 10 ** (dB_values / 10)
#noise_levels = np.array([1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20])
#noise_levels = 10**(noise_levels/10)
#noise_levels = [0.001, 0.002, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
#noise_levels = [1000, 500, 100, 50, 20, 10, 5, 2, 1]
#noise_levels = [10, 20, 50] #the previous correspondence for this was 0.0007
d_bers = []
for noise_level in noise_levels:
    d_ber = transmit(5, noise_level)
    print(N, d_ber)
    d_bers.append(d_ber)
print(d_bers)
Eb_N0_dB = noise_levels
theoryBerAWGN = 0.5 * erfc(np.sqrt(10**(Eb_N0_dB / 10)))
EbN0Lin = 10**(Eb_N0_dB / 10)
theoryBerRayleigh = 0.5 * (1 - np.sqrt(EbN0Lin / (EbN0Lin + 1)))

plt.figure()
#plt.semilogy(Eb_N0_dB, theoryBerAWGN, 'cd-', linewidth=2, label='AWGN Theory')
plt.semilogy(Eb_N0_dB, theoryBerRayleigh, 'bp-', linewidth=2, label='Rayleigh Theory')
plt.semilogy(noise_levels, d_bers, 'mx-',label = 'Data BERs', marker = 'o')
plt.legend()
plt.xlabel("Noise level in (1/SNR) in linear scale")
plt.ylabel('Error value')
def log10_tick_formatter(val, pos=None):
    return f'{10 * np.log10(val):.1f}'
#plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(log10_tick_formatter))
plt.title("Plot of MSE, BER for fixed doppler frequency 100Hz as a function of (1/SNR)")
plt.show()
