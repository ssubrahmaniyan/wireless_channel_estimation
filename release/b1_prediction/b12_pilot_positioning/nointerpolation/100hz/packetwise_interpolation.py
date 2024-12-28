import numpy as np
import json
import matplotlib.pyplot as plt

# Load the complex channel values
with open("channels.txt", "r") as f:
    channel_vals = [json.loads(line)[0] for line in f]
channel_vals = np.array([complex(re, im) for re, im in channel_vals])

def generate_random_bits(n):
    return np.random.choice([-1, 1], size=n)

def transmit_data(data, channel, snr_db):
    signal = channel * data
    signal_power = np.mean(np.abs(signal)**2)
    noise_power = signal_power / (10**(snr_db/10))
    noise = np.sqrt(noise_power/2) * (np.random.randn(len(data)) + 1j * np.random.randn(len(data)))
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

def interpolate_channel(pilot_indices, channel_estimates, current_index):
    lower_index = max([idx for idx in pilot_indices if idx <= current_index])
    upper_index = min([idx for idx in pilot_indices if idx > current_index], default=lower_index)

    if lower_index == upper_index:
        return channel_estimates[pilot_indices.index(lower_index)]

    h_lower = channel_estimates[pilot_indices.index(lower_index)]
    h_upper = channel_estimates[pilot_indices.index(upper_index)]
    alpha = (current_index - lower_index) / (upper_index - lower_index)
    return (1 - alpha) * h_lower + alpha * h_upper

def run_simulation(snr_db, mse_threshold, packet_size=50, pilot_size=50):
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
    
    channel_estimates.append(h_est)
    pilot_indices.append(current_index - pilot_size)
    total_bits += pilot_size

    # Main loop
    while current_index < len(channel_vals) - packet_size:
        # Interpolate channel
        h_interp = interpolate_channel(pilot_indices, channel_estimates, current_index)

        # Data transmission
        data_bits = generate_random_bits(packet_size)
        actual_channel = channel_vals[current_index:current_index + packet_size]
        received_data = transmit_data(data_bits, actual_channel, snr_db)

        # Demodulate using interpolated channel estimate
        demodulated_bits = np.sign(np.real(received_data * np.conjugate(h_interp) / (np.abs(h_interp)**2)))

        # Calculate bit errors
        errors = np.sum(np.abs(demodulated_bits - data_bits) / 2)
        bit_errors += errors
        total_data_bits += packet_size

        # Re-estimate channel
        h_new = estimate_channel(received_data, demodulated_bits)

        # Calculate MSE
        mse = np.abs(h_new - h_interp)**2
        mse_values.append(mse)

        if mse > mse_threshold:
            # Send new pilot
            h_est, current_index = send_pilot(current_index, pilot_size, snr_db)
            if h_est is None:
                break
            channel_estimates.append(h_est)
            pilot_indices.append(current_index - pilot_size)
            total_bits += pilot_size
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
snr_values = np.arange(0, 21, 2)  # 0 to 20 dB in steps of 2
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
            bers.append(ber)
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

# Create plots
plt.figure(figsize=(15, 6))
colors = ['b', 'g', 'r', 'c', 'm']
markers = ['o', 's', '^', 'D', 'v']

# Pilot Ratio Plot
plt.subplot(1, 2, 1)
for (threshold, ratios), color, marker in zip(results.items(), colors, markers):
    plt.plot(snr_values, ratios, f'{color}-{marker}', label=f'MSE Threshold = {threshold}')

plt.xlabel('SNR (dB)')
plt.ylabel('Pilot Ratio')
plt.title('Pilot Ratio vs SNR (Interpolation)')
plt.grid(True)
plt.legend()

# BER Plot
plt.subplot(1, 2, 2)
for threshold, color, marker in zip(mse_thresholds, colors, markers):
    bers = all_metrics[threshold]['bers']
    plt.plot(snr_values, bers, f'{color}-{marker}', label=f'MSE Threshold = {threshold}')

plt.xlabel('SNR (dB)')
plt.ylabel('Bit Error Rate')
plt.title('BER vs SNR (Interpolation)')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.savefig('simulation_results_interpolation.png', dpi=300, bbox_inches='tight')
plt.show()

# Save detailed results to text file
with open('simulation_results_interpolation.txt', 'w') as f:
    f.write("Channel Estimation Simulation Results (Interpolation)\n")
    f.write("==============================================\n\n")
    
    for threshold in mse_thresholds:
        f.write(f"\nMSE Threshold = {threshold}\n")
        f.write("SNR (dB) | Pilot Ratio |    BER     | Avg MSE\n")
        f.write("-" * 50 + "\n")
        
        for i, snr in enumerate(snr_values):
            metrics = all_metrics[threshold]
            f.write(f"{snr:8.1f} | {metrics['pilot_ratios'][i]:10.4f} | {metrics['bers'][i]:10.4f} | {metrics['avg_mses'][i]:8.4f}\n")
        
        f.write("\n" + "=" * 50 + "\n")

# Save as CSV for easier data processing
with open('simulation_results_interpolation.csv', 'w') as f:
    # Header
    f.write("mse_threshold,snr,pilot_ratio,ber,avg_mse\n")
    
    # Data
    for threshold in mse_thresholds:
        for i, snr in enumerate(snr_values):
            metrics = all_metrics[threshold]
            f.write(f"{threshold},{snr},{metrics['pilot_ratios'][i]},{metrics['bers'][i]},{metrics['avg_mses'][i]}\n")

print("Results have been saved to:")
print("1. simulation_results_interpolation.txt - Formatted text file")
print("2. simulation_results_interpolation.csv - CSV file for data processing")
print("3. simulation_results_interpolation.png - Plot image with both Pilot Ratio and BER")
