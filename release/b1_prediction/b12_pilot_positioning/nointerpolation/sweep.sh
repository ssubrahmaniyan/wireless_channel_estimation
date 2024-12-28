#!/bin/bash

# Create or clear the output files
echo "dopp,snr,ratio" > snr_sweep_data.txt
echo "dopp,threshold,ratio" > mse_sweep_data.txt

# Define parameter ranges
doppler_freqs=($(seq 800 -60 10))  # Doppler frequencies from 800 to 10
snr_values=(30 20 10 5 1)   # SNR values
mse_thresholds=(0.01 0.05 0.1 0.2 0.5)  # MSE threshold values

# Fixed parameters
fixed_snr=10
fixed_mse=0.1

# First sweep: Different SNRs at fixed MSE threshold
echo "Running SNR sweep..."
for dopp in "${doppler_freqs[@]}"; do
    echo "Processing Doppler frequency: $dopp Hz"
    # Generate channel values
    python3 jakeschannel.py -T comp -Fd "$dopp" -P 100000
    
    # Run simulations for each SNR
    for snr in "${snr_values[@]}"; do
        # Run simulation and extract pilot ratio
        result=$(python3 - <<EOF
import numpy as np
import json

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

    while current_index < len(channel_vals) - packet_size:
        data_bits = generate_random_bits(packet_size)
        actual_channel = channel_vals[current_index:current_index + packet_size]
        received_data = transmit_data(data_bits, actual_channel, snr_db)
        
        h_prev = channel_estimates[-1]
        demodulated_bits = np.sign(np.real(received_data * np.conjugate(h_prev) / (np.abs(h_prev)**2)))
        
        errors = np.sum(np.abs(demodulated_bits - data_bits) / 2)
        bit_errors += errors
        total_data_bits += packet_size
        
        h_new = estimate_channel(received_data, demodulated_bits)
        mse = np.abs(h_new - h_prev)**2
        mse_values.append(mse)
        
        if mse > mse_threshold:
            h_est, current_index = send_pilot(current_index, pilot_size, snr_db)
            if h_est is None:
                break
            channel_estimates.append(h_est)
            pilot_indices.append(current_index - pilot_size)
            total_bits += pilot_size
        else:
            channel_estimates.append(h_new)
            current_index += packet_size
            total_bits += packet_size

    pilot_ratio = (len(pilot_indices) * pilot_size) / total_bits if total_bits > 0 else 0
    return pilot_ratio

ratio = run_simulation($snr, $fixed_mse)
print(ratio)
EOF
)
        echo "$dopp,$snr,$result" >> snr_sweep_data.txt
    done
done

# Second sweep: Different MSE thresholds at fixed SNR
echo "Running MSE threshold sweep..."
for dopp in "${doppler_freqs[@]}"; do
    echo "Processing Doppler frequency: $dopp Hz"
    python3 jakeschannel.py -T comp -Fd "$dopp" -P 300000
    
    for threshold in "${mse_thresholds[@]}"; do
        result=$(python3 - <<EOF
import numpy as np
import json

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

    while current_index < len(channel_vals) - packet_size:
        data_bits = generate_random_bits(packet_size)
        actual_channel = channel_vals[current_index:current_index + packet_size]
        received_data = transmit_data(data_bits, actual_channel, snr_db)
        
        h_prev = channel_estimates[-1]
        demodulated_bits = np.sign(np.real(received_data * np.conjugate(h_prev) / (np.abs(h_prev)**2)))
        
        errors = np.sum(np.abs(demodulated_bits - data_bits) / 2)
        bit_errors += errors
        total_data_bits += packet_size
        
        h_new = estimate_channel(received_data, demodulated_bits)
        mse = np.abs(h_new - h_prev)**2
        mse_values.append(mse)
        
        if mse > mse_threshold:
            h_est, current_index = send_pilot(current_index, pilot_size, snr_db)
            if h_est is None:
                break
            channel_estimates.append(h_est)
            pilot_indices.append(current_index - pilot_size)
            total_bits += pilot_size
        else:
            channel_estimates.append(h_new)
            current_index += packet_size
            total_bits += packet_size

    pilot_ratio = (len(pilot_indices) * pilot_size) / total_bits if total_bits > 0 else 0
    return pilot_ratio
ratio = run_simulation($fixed_snr, $threshold)
print(ratio)
EOF
)
        echo "$dopp,$threshold,$result" >> mse_sweep_data.txt
    done
done

# Create plots
python3 - <<EOF
import matplotlib.pyplot as plt
import pandas as pd

# Plot 1: SNR sweep
df_snr = pd.read_csv('snr_sweep_data.txt')
plt.figure(figsize=(12, 6))

for snr in sorted(df_snr['snr'].unique(), reverse=True):
    data = df_snr[df_snr['snr'] == snr]
    plt.plot(data['dopp'], data['ratio'], marker='o', label=f'SNR = {snr} dB', markersize=4)

plt.xlabel('Doppler Frequency (Hz)')
plt.ylabel('Pilot Ratio')
plt.title('Pilot Ratio vs Doppler Frequency for Different SNRs (MSE = 0.1)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('snr_sweep.png', dpi=300, bbox_inches='tight')
plt.close()

# Plot 2: MSE threshold sweep
df_mse = pd.read_csv('mse_sweep_data.txt')
plt.figure(figsize=(12, 6))

for threshold in sorted(df_mse['threshold'].unique()):
    data = df_mse[df_mse['threshold'] == threshold]
    plt.plot(data['dopp'], data['ratio'], marker='o', label=f'MSE = {threshold}', markersize=4)

plt.xlabel('Doppler Frequency (Hz)')
plt.ylabel('Pilot Ratio')
plt.title('Pilot Ratio vs Doppler Frequency for Different MSE Thresholds (SNR = 10 dB)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('mse_sweep.png', dpi=300, bbox_inches='tight')
plt.close()
EOF

echo "Simulation complete. Results saved as:"
echo "1. snr_sweep.png - Plot for different SNR values"
echo "2. mse_sweep.png - Plot for different MSE thresholds"
echo "Raw data saved in snr_sweep_data.txt and mse_sweep_data.txt"