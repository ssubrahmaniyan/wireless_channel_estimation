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
    signal = channel / np.sqrt(channel_mod) * data
    eb_no_linear = 10**(eb_no_db / 10)
    n0 = 1 / eb_no_linear
    noise = np.sqrt(n0 / 2) * (np.random.randn(len(data)) + 1j * np.random.randn(len(data)))
    return signal + noise

def estimate_channel(received, data):
    X = data.reshape(-1, 1)
    h = np.linalg.lstsq(X, received, rcond=None)[0][0]
    return h

def process_pilot(current_index, pilot_size, snr_db):
    if current_index + pilot_size > len(channel_vals):
        return None, None, None, current_index
    
    pilot_bits = generate_random_bits(pilot_size)
    actual_channel = channel_vals[current_index:current_index + pilot_size]
    received_pilot = transmit_data(pilot_bits, actual_channel, snr_db)
    h_est = estimate_channel(received_pilot, pilot_bits)
    
    return h_est, pilot_bits, received_pilot, current_index + pilot_size

def linear_interpolation(start_val, end_val, num_points):
    if num_points <= 1:
        return np.array([start_val])
    t = np.linspace(0, 1, num_points)
    return (1 - t[:, np.newaxis]) * start_val + t[:, np.newaxis] * end_val

def run_simulation(snr_db, mse_threshold, packet_size=5, pilot_size=5):
    current_index = 0
    total_bits = 0
    bit_errors = 0
    total_data_bits = 0
    data_buffer = []  # Store data between pilots for reprocessing
    
    # Initial pilot
    h_pilot_start, pilot_bits, received_pilot, current_index = process_pilot(current_index, pilot_size, snr_db)
    if h_pilot_start is None:
        return None
    
    total_bits += pilot_size
    h_current = h_pilot_start
    consecutive_packets = 0
    
    while current_index < len(channel_vals) - packet_size:
        # Data transmission
        data_bits = generate_random_bits(packet_size)
        actual_channel = channel_vals[current_index:current_index + packet_size]
        received_data = transmit_data(data_bits, actual_channel, snr_db)
        
        # Demodulate using current channel estimate
        demodulated = received_data * np.conjugate(h_current) / (np.abs(h_current)**2)
        demodulated_bits = np.array([1 if np.real(x) > 0 else -1 for x in demodulated])
        
        # Store data for possible reprocessing
        data_buffer.append((data_bits, received_data, current_index))
        consecutive_packets += 1
        
        # Re-estimate channel using demodulated data
        h_new = estimate_channel(received_data, demodulated_bits)
        
        # Calculate MSE
        mse = np.abs(h_new - h_current)**2
        
        if mse > mse_threshold or consecutive_packets >= 5:
            # Send new pilot
            h_pilot_end, new_pilot_bits, new_received_pilot, next_index = process_pilot(
                current_index + packet_size, pilot_size, snr_db)
            
            if h_pilot_end is None:
                break
                
            # Interpolate channel estimates
            num_packets = len(data_buffer)
            interpolated_channels = linear_interpolation(h_pilot_start, h_pilot_end, num_packets)
            
            # Reprocess stored data with interpolated channels
            for idx, (stored_bits, stored_received, stored_index) in enumerate(data_buffer):
                h_interp = interpolated_channels[idx]
                demod_with_interp = stored_received * np.conjugate(h_interp) / (np.abs(h_interp)**2)
                final_bits = np.array([1 if np.real(x) > 0 else -1 for x in demod_with_interp])
                bit_errors += np.sum(np.abs(final_bits - stored_bits) / 2)
                total_data_bits += len(stored_bits)
            
            # Reset for next segment
            h_pilot_start = h_pilot_end
            current_index = next_index
            data_buffer = []
            consecutive_packets = 0
            total_bits += pilot_size
            
        else:
            h_current = h_new
            current_index += packet_size
            total_bits += packet_size
    
    # Process any remaining data in buffer
    if data_buffer:
        for stored_bits, stored_received, stored_index in data_buffer:
            demod = stored_received * np.conjugate(h_current) / (np.abs(h_current)**2)
            final_bits = np.array([1 if np.real(x) > 0 else -1 for x in demod])
            bit_errors += np.sum(np.abs(final_bits - stored_bits) / 2)
            total_data_bits += len(stored_bits)
    
    # Calculate final metrics
    ber = bit_errors / total_data_bits if total_data_bits > 0 else 0
    pilot_ratio = (total_bits - total_data_bits) / total_bits
    
    return pilot_ratio, ber

# Parameter sweep
snr_values = np.arange(0, 21, 2)
mse_thresholds = [0.01, 0.05, 0.1, 0.2, 0.5]
results = {}

# Run simulations
for threshold in mse_thresholds:
    pilot_ratios = []
    bers = []
    for snr in snr_values:
        result = run_simulation(snr, threshold)
        if result is not None:
            pilot_ratio, ber = result
            pilot_ratios.append(pilot_ratio)
            bers.append(ber)
        else:
            pilot_ratios.append(np.nan)
            bers.append(np.nan)
    results[threshold] = {
        'pilot_ratios': pilot_ratios,
        'bers': bers
    }

# Add theoretical curves
def q_function(x):
    return 0.5 * erfc(x/np.sqrt(2))

def theoretical_bpsk_awgn(snr_db):
    snr_linear = 10**(snr_db/10)
    return q_function(np.sqrt(2*snr_linear))

def theoretical_bpsk_perfect_csi(snr_db):
    # For perfect CSI, performance is same as AWGN but with 3dB loss due to normalization
    return theoretical_bpsk_awgn(snr_db - 3)

def theoretical_bpsk_imperfect_csi(snr_db, pilot_len=5):
    # Model imperfect CSI with estimation error
    snr_linear = 10**(snr_db/10)
    pilot_snr = snr_linear * pilot_len  # SNR improvement due to pilot length
    variance_h_est = 1/(1 + pilot_snr)  # Channel estimation error variance
    effective_snr = snr_linear / (1 + variance_h_est * snr_linear)
    return theoretical_bpsk_awgn(10*np.log10(effective_snr) - 3)

def theoretical_bpsk_rayleigh(snr_db):
    # For Rayleigh fading with perfect CSI
    snr_linear = 10**(snr_db/10)
    return 0.5 * (1 - np.sqrt(snr_linear/(1 + snr_linear)))

# [Previous theoretical functions remain the same]

# Create plots with theoretical curves including Rayleigh
plt.figure(figsize=(15, 6))
colors = ['b', 'g', 'r', 'c', 'm']
markers = ['o', 's', '^', 'D', 'v']

# Pilot Ratio Plot [remains the same]
plt.subplot(1, 2, 1)
for (threshold, metrics), color, marker in zip(results.items(), colors, markers):
    plt.plot(snr_values, metrics['pilot_ratios'], f'{color}-{marker}', 
             label=f'MSE Threshold = {threshold}')

plt.xlabel('SNR (dB)')
plt.ylabel('Pilot Ratio')
plt.title('Pilot Ratio vs SNR (Interpolation)')
plt.grid(True)
plt.legend()

# BER Plot with theoretical curves including Rayleigh
plt.subplot(1, 2, 2)
# Plot simulation results
for (threshold, metrics), color, marker in zip(results.items(), colors, markers):
    plt.semilogy(snr_values, metrics['bers'], f'{color}-{marker}', 
                 label=f'MSE Threshold = {threshold}')

# Add theoretical curves
theoretical_snr = np.linspace(0, 20, 100)
plt.semilogy(theoretical_snr, 
             [theoretical_bpsk_awgn(snr) for snr in theoretical_snr],
             'k--', label='BPSK AWGN (Theory)')
plt.semilogy(theoretical_snr,
             [theoretical_bpsk_perfect_csi(snr) for snr in theoretical_snr],
             'k:', label='BPSK Perfect CSI (Theory)')
plt.semilogy(theoretical_snr,
             [theoretical_bpsk_imperfect_csi(snr) for snr in theoretical_snr],
             'k-.', label='BPSK Imperfect CSI (Theory)')
plt.semilogy(theoretical_snr,
             [theoretical_bpsk_rayleigh(snr) for snr in theoretical_snr],
             'r--', label='BPSK Rayleigh (Theory)')

plt.xlabel('SNR (dB)')
plt.ylabel('Bit Error Rate')
plt.title('BER vs SNR (Interpolation)')
plt.grid(True)
plt.legend()
plt.ylim([1e-5, 1])

plt.tight_layout()
plt.savefig('simulation_results_interpolation_with_theory.png', dpi=300, bbox_inches='tight')
plt.show()

# Save results with theoretical values including Rayleigh
with open('simulation_results_interpolation_with_theory.csv', 'w') as f:
    f.write("snr,awgn_theory,perfect_csi_theory,imperfect_csi_theory,rayleigh_theory\n")
    for snr in theoretical_snr:
        awgn = theoretical_bpsk_awgn(snr)
        perfect = theoretical_bpsk_perfect_csi(snr)
        imperfect = theoretical_bpsk_imperfect_csi(snr)
        rayleigh = theoretical_bpsk_rayleigh(snr)
        f.write(f"{snr},{awgn},{perfect},{imperfect},{rayleigh}\n")