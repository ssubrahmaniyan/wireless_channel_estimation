import numpy as np
import json
import matplotlib.pyplot as plt
from scipy.special import erfc 

# Load the complex channel values
with open("channels.txt", "r") as f:
    channel_vals = [json.loads(line)[0] for line in f]
channel_vals = np.array([complex(re, im) for re, im in channel_vals])
channel_mod = np.mean(np.abs(channel_vals) ** 2)

def generate_random_bits(n):
    return np.random.choice([-1, 1], size=n)

def transmit_data(data, channel, eb_no_db, bit_rate=1):
    #signal = channel * data
    signal = channel / np.sqrt(channel_mod) * data
    eb_no_linear = 10**(eb_no_db / 10)
    n0 = 1 / eb_no_linear  # Fixed noise power calculation
    noise = np.sqrt(n0 / 2) * (np.random.randn(len(data)) + 1j * np.random.randn(len(data)))
    return signal + noise

def estimate_channel(received, data):
    X = data.reshape(-1, 1)
    h = np.linalg.lstsq(X, received, rcond=None)[0][0]
    return h

def send_pilot(current_index, pilot_size, snr_db):
    if current_index + pilot_size > len(channel_vals):
        return None, current_index
        
    pilot_bits = generate_random_bits(pilot_size)
    actual_channel = channel_vals[current_index:current_index + pilot_size]
    received_pilot = transmit_data(pilot_bits, actual_channel, snr_db)
    h_est = estimate_channel(received_pilot, pilot_bits)
    return h_est, current_index + pilot_size

def sign(number):
    if np.angle(number) > -np.pi/2 and np.angle(number) < np.pi/2:
        return 1
    return -1

def run_simulation(snr_db, mse_threshold, packet_size=5, pilot_size=5):
    current_index = 0
    channel_estimates = []
    pilot_indices = []
    total_bits = 0
    bit_errors = 0
    total_data_bits = 0
    mse_values = []
    index = 0
    # Initial pilot
    h_est, current_index = send_pilot(current_index, pilot_size, snr_db) # Sending y = hp + n
    if h_est is None:
        return None
    
    channel_estimates.append(h_est) # Finding h from y = hp + n and appending it
    pilot_indices.append(current_index - pilot_size)
    total_bits += pilot_size
    currmse = []
    # Main loop
    while current_index < len(channel_vals) - packet_size:
        index += 1
        # Data transmission
        data_bits = generate_random_bits(packet_size) # generate data packet [x]
        actual_channel = channel_vals[current_index:current_index + packet_size] # find the vector array [h']
        received_data = transmit_data(data_bits, actual_channel, snr_db) # return [y] = [h'][x] + [n]
        
        # Demodulate using previous channel estimate
        h_prev = channel_estimates[-1]
        demodulated = received_data * np.conjugate(h_prev) / (np.abs(h_prev)**2)
        demodulated_bits = np.array([1 if np.real(x) > 0 else -1 for x in demodulated])
        
        # Re-estimate channel
        h_new = np.linalg.lstsq(data_bits.reshape(-1, 1), received_data, rcond=None)[0][0]        
        channel_estimates.append(h_new)

        # Re-do the demodulation of data
        demodulated = received_data * np.conjugate(h_new) / (np.abs(h_new)**2)
        demodulated_bits = np.array([1 if np.real(x) > 0 else -1 for x in demodulated])        # plt.plot(data_bits, label = "data")
        # plt.plot(demodulated_bits, label = "demod")
        # plt.plot(actual_channel.real,label = "real channel")
        # plt.plot(actual_channel.imag, label = 'imag channel')
        # plt.legend()
        # plt.show()
        errors = 0
        for i in range(len(demodulated_bits)):
            if demodulated_bits[i] != data_bits[i]:
                errors += 1
        bit_errors += errors
        total_data_bits += packet_size
        total_bits += packet_size
        
        # print(sign(h_prev), sign(h_new))
        # plt.plot(data_bits, label = 'original')
        # plt.plot(demodulated_bits, label='demodulated bits')
        # plt.show()

        # Calculate MSE
        mse = np.abs(h_new - h_prev)**2
        mse_values.append(mse)
        # currmse.append(mse)
        # avgmse = np.mean(currmse)
        if mse > mse_threshold:
        # if index == 1:
            # Send new pilot
            h_est, current_index = send_pilot(current_index, pilot_size, snr_db)
            if h_est is None:
                break
            channel_estimates.append(h_est)
            pilot_indices.append(current_index - pilot_size)
            total_bits += pilot_size
            currmse = []
            index = 0
        else:
            # Continue with data
            current_index += packet_size
            total_bits += packet_size

    # Calculate metrics
    ber = bit_errors / total_data_bits if total_data_bits > 0 else 0
    pilot_ratio = (len(pilot_indices) * pilot_size) / total_bits
    avg_mse = np.mean(mse_values) if mse_values else 0
    
    return pilot_ratio, ber, avg_mse

# Parameter sweep
snr_values = np.arange(0,40, 5)  # 0 to 30 dB in steps of 5 dB each
mse_thresholds = [0.01, 0.05, 0.1, 0.2, 0.5]
results = {}
all_metrics = {}

# Run simulations
for threshold in mse_thresholds:
    pilot_ratios = []
    bers = []
    avg_mses = []
    for snr in snr_values:
        result = run_simulation(snr, threshold)
        if result is not None:
            pilot_ratio, ber, avg_mse = result
            pilot_ratios.append(pilot_ratio)
            if ber != 0:
                bers.append(ber)
            else:
                bers.append(1e-6)
            avg_mses.append(avg_mse)
        else:
            pilot_ratios.append(np.nan)
            bers.append(np.nan)
            avg_mses.append(np.nan)
    results[threshold] = pilot_ratios
    all_metrics[threshold] = {
        'pilot_ratios': pilot_ratios,
        'bers': bers,
        'avg_mses': avg_mses
    }

plt.figure(figsize=(20, 6))
colors = ['b', 'g', 'r', 'c', 'm']
markers = ['o', 's', '^', 'D', 'v']

# Pilot Ratio Plot
plt.subplot(1, 3, 1)
for (threshold, ratios), color, marker in zip(results.items(), colors, markers):
    plt.plot(snr_values, ratios, f'{color}-{marker}', label=f'MSE Threshold = {threshold}')
plt.xlabel('SNR (dB)')
plt.ylabel('Pilot Ratio')
plt.title('Pilot Ratio vs SNR')
plt.grid(True)
plt.legend()

# BER Plot
plt.subplot(1, 3, 2)
# Calculate theoretical BER curves
EbN0Lin = 10**(snr_values / 10)
theoryBerRayleigh = 0.5 * (1 - np.sqrt(EbN0Lin / (EbN0Lin + 1)))

# Plot theoretical curves
plt.semilogy(snr_values, theoryBerRayleigh, 'k:', linewidth=2, label='Rayleigh Theory')

# Plot simulation results
for threshold, color, marker in zip(mse_thresholds, colors, markers):
    bers = all_metrics[threshold]['bers']
    plt.semilogy(snr_values, bers, f'{color}-{marker}', label=f'MSE Threshold = {threshold}')

plt.xlabel('SNR (dB)')
plt.ylabel('Bit Error Rate (BER)')
plt.title('BER vs SNR')
plt.grid(True)
plt.legend()
plt.ylim([1e-5, 0.5])

# MSE Plot
plt.subplot(1, 3, 3)
for threshold, color, marker in zip(mse_thresholds, colors, markers):
    avg_mses = all_metrics[threshold]['avg_mses']
    plt.plot(snr_values, avg_mses, f'{color}-{marker}', label=f'MSE Threshold = {threshold}')

plt.xlabel('SNR (dB)')
plt.ylabel('Average MSE')
plt.title('Average MSE vs SNR')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.savefig('simulation_results_with_mse.png', dpi=300, bbox_inches='tight')
plt.show()

# Save detailed results to text file
with open('simulation_results.txt', 'w') as f:
    f.write("Channel Estimation Simulation Results\n")
    f.write("===================================\n\n")
    
    for threshold in mse_thresholds:
        f.write(f"\nMSE Threshold = {threshold}\n")
        f.write("SNR (dB) | Pilot Ratio |    BER     | Avg MSE\n")
        f.write("-" * 50 + "\n")
        
        for i, snr in enumerate(snr_values):
            metrics = all_metrics[threshold]
            f.write(f"{snr:8.1f} | {metrics['pilot_ratios'][i]:10.4f} | {metrics['bers'][i]:10.4f} | {metrics['avg_mses'][i]:8.4f}\n")
        
        f.write("\n" + "=" * 50 + "\n")

# Save as CSV for easier data processing
with open('simulation_results.csv', 'w') as f:
    # Header
    f.write("mse_threshold,snr,pilot_ratio,ber,avg_mse\n")
    
    # Data
    for threshold in mse_thresholds:
        for i, snr in enumerate(snr_values):
            metrics = all_metrics[threshold]
            f.write(f"{threshold},{snr},{metrics['pilot_ratios'][i]},{metrics['bers'][i]},{metrics['avg_mses'][i]}\n")

print("Results have been saved to:")
print("1. simulation_results.txt - Formatted text file")
print("2. simulation_results.csv - CSV file for data processing")
print("3. simulation_results.png - Plot image with both Pilot Ratio and BER")