import numpy as np
import json
import matplotlib.pyplot as plt
from scipy.special import erfc
from scipy.interpolate import interp1d

# Load the complex channel values
with open("channels.txt", "r") as f:
    channel_vals = [json.loads(line)[0] for line in f]
channel_vals = np.array([complex(re, im) for re, im in channel_vals])
channel_mod = np.mean(np.abs(channel_vals) ** 2)

def generate_random_bits(n):
    return np.random.choice([-1, 1], size=n)

def transmit_data(data, channel, eb_no_db, bit_rate=1):
    signal = channel / np.sqrt(channel_mod) * data
    eb_no_linear = 10**(eb_no_db / 10)
    n0 = 1 / eb_no_linear
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

def interpolate_channel(time_points, channel_estimates):
    # Interpolate real and imaginary parts separately
    real_interp = interp1d(time_points, np.real(channel_estimates), kind='linear')
    imag_interp = interp1d(time_points, np.imag(channel_estimates), kind='linear')
    
    # Return a function that combines both interpolations
    return lambda x: real_interp(x) + 1j * imag_interp(x)

def run_simulation(snr_db, mse_threshold, packet_size=5, pilot_size=5):
    current_index = 0
    channel_estimates = []
    pilot_indices = []
    total_bits = 0
    bit_errors = 0
    total_data_bits = 0
    mse_values = []
    
    # Initial pilot
    h_est, current_index = send_pilot(current_index, pilot_size, snr_db)
    if h_est is None:
        return None
    
    channel_0 = h_est
    pilot_time = current_index - pilot_size + pilot_size//2  # Center of pilot
    
    while current_index < len(channel_vals) - 2*packet_size:  # Need space for two packets
        # First data packet
        data_bits_1 = generate_random_bits(packet_size)
        actual_channel_1 = channel_vals[current_index:current_index + packet_size]
        received_data_1 = transmit_data(data_bits_1, actual_channel_1, snr_db)
        
        # Demodulate first packet using pilot channel estimate
        demodulated_1 = received_data_1 * np.conjugate(channel_0) / (np.abs(channel_0)**2)
        data_1 = np.array([1 if np.real(x) > 0 else -1 for x in demodulated_1])
        
        # Estimate channel for first packet
        channel_1A = estimate_channel(received_data_1, data_1)
        time_1 = current_index + packet_size//2
        
        # Second data packet
        current_index += packet_size
        data_bits_2 = generate_random_bits(packet_size)
        actual_channel_2 = channel_vals[current_index:current_index + packet_size]
        received_data_2 = transmit_data(data_bits_2, actual_channel_2, snr_db)
        
        # Demodulate second packet using channel_1A
        demodulated_2 = received_data_2 * np.conjugate(channel_1A) / (np.abs(channel_1A)**2)
        data_2 = np.array([1 if np.real(x) > 0 else -1 for x in demodulated_2])
        
        # Estimate channel for second packet
        channel_2A = estimate_channel(received_data_2, data_2)
        time_2 = current_index + packet_size//2
        
        # Interpolate channel for first packet
        interp_func = interpolate_channel(
            [pilot_time, time_2],
            [channel_0, channel_2A]
        )
        interpolated_channel = interp_func(time_1)
        
        # Re-demodulate first packet using interpolated channel
        demodulated_1_final = received_data_1 * np.conjugate(interpolated_channel) / (np.abs(interpolated_channel)**2)
        final_bits_1 = np.array([1 if np.real(x) > 0 else -1 for x in demodulated_1_final])
        
        # Calculate errors for first packet
        errors = np.sum(final_bits_1 != data_bits_1)
        bit_errors += errors
        total_data_bits += packet_size
        
        # Calculate MSE between successive channel estimates
        mse = np.abs(channel_2A - channel_0)**2
        mse_values.append(mse)
        
        if mse > mse_threshold:
            # Send new pilot
            h_est, current_index = send_pilot(current_index, pilot_size, snr_db)
            if h_est is None:
                break
            channel_0 = h_est
            pilot_time = current_index - pilot_size + pilot_size//2
            pilot_indices.append(current_index - pilot_size)
            total_bits += pilot_size
        else:
            # Continue with next iteration
            channel_0 = channel_2A
            pilot_time = time_2
            current_index += packet_size
            
        total_bits += 2 * packet_size  # Count both packets
    
    # Calculate metrics
    ber = bit_errors / total_data_bits if total_data_bits > 0 else 0
    pilot_ratio = (len(pilot_indices) * pilot_size) / total_bits
    avg_mse = np.mean(mse_values) if mse_values else 0
    
    return pilot_ratio, ber, avg_mse

# Parameter sweep
snr_values = np.arange(0, 40, 5)
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

# Plotting code remains the same
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

plt.semilogy(snr_values, theoryBerRayleigh, 'k:', linewidth=2, label='Rayleigh Theory')

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
plt.savefig('simulation_results_with_interpolation.png', dpi=300, bbox_inches='tight')
plt.show()

# Save results to files
with open('simulation_results_interpolation.txt', 'w') as f:
    f.write("Channel Estimation Simulation Results (with Interpolation)\n")
    f.write("================================================\n\n")
    
    for threshold in mse_thresholds:
        f.write(f"\nMSE Threshold = {threshold}\n")
        f.write("SNR (dB) | Pilot Ratio |    BER     | Avg MSE\n")
        f.write("-" * 50 + "\n")
        
        for i, snr in enumerate(snr_values):
            metrics = all_metrics[threshold]
            f.write(f"{snr:8.1f} | {metrics['pilot_ratios'][i]:10.4f} | {metrics['bers'][i]:10.4f} | {metrics['avg_mses'][i]:8.4f}\n")
        
        f.write("\n" + "=" * 50 + "\n")

with open('simulation_results_interpolation.csv', 'w') as f:
    f.write("mse_threshold,snr,pilot_ratio,ber,avg_mse\n")
    
    for threshold in mse_thresholds:
        for i, snr in enumerate(snr_values):
            metrics = all_metrics[threshold]
            f.write(f"{threshold},{snr},{metrics['pilot_ratios'][i]},{metrics['bers'][i]},{metrics['avg_mses'][i]}\n")

print("Results have been saved to:")
print("1. simulation_results_interpolation.txt - Formatted text file")
print("2. simulation_results_interpolation.csv - CSV file for data processing")
print("3. simulation_results_with_interpolation.png - Plot image")