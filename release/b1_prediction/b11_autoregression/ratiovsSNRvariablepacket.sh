#!/bin/bash

# Create or clear the output file
output_file="ratio_vs_SNR_data.txt"

echo "dopp,SNR,ratio" > $output_file

# Define parameter ranges
doppler_freqs=($(seq 500 -100 10))  # Doppler frequencies from 1000 to 100
snr_thresholds=(30 20 10 5 1)   # MSE thresholds

# Loop over all parameters
for snr in "${snr_thresholds[@]}"
do
    for dopp in "${doppler_freqs[@]}"
    do
        # Run Jake's channel model
        python3 jakeschannel.py -T comp -Fd "$dopp" -P 500000
        
        # Run the channel estimation with current parameters and capture the ratio
        ratio=$(python3 - <<END
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.tsa.api import VAR
import json

# Load the complex channel values
with open("channels.txt", "r") as f:
    channel_vals = [json.loads(line)[0] for line in f]
channel_vals = np.array([complex(re, im) for re, im in channel_vals])

# Parameters
packet_size = int((1/20) * (100000/$dopp))
pilot_size = packet_size
var_order = int(packet_size/2)
mse_threshold = 0.1
snr_db = $snr
initial_pilots = var_order + 2

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

def collect_pilot_estimates(current_index, num_pilots):
    estimates = []
    for _ in range(num_pilots):
        if current_index + pilot_size > len(channel_vals):
            break
        pilot_bits = generate_random_bits(pilot_size)
        actual_channel = channel_vals[current_index:current_index + pilot_size]
        received_pilot = transmit_data(pilot_bits, actual_channel, snr_db)
        hbar = estimate_channel(received_pilot, pilot_bits)
        estimates.append(hbar)
        current_index += pilot_size
    return estimates, current_index

def predict_channel_var(channel_history, steps=1):
    if len(channel_history) <= var_order:
        raise ValueError(f"Not enough samples for VAR")
    channel_history_separated = np.column_stack((np.real(channel_history), np.imag(channel_history)))
    model = VAR(channel_history_separated)
    results = model.fit(var_order)
    forecast = results.forecast(channel_history_separated[-var_order:], steps=steps)
    return forecast[0, 0] + 1j * forecast[0, 1]

current_index = 0
channel_estimates = []
pilot_requests = []
total_bits = 0

initial_estimates, current_index = collect_pilot_estimates(current_index, initial_pilots)
channel_estimates.extend(initial_estimates)
pilot_requests.extend(range(0, current_index, pilot_size))
total_bits += initial_pilots * pilot_size

while current_index < len(channel_vals) - packet_size:
    try:
        htilde = predict_channel_var(channel_estimates)
    except ValueError:
        break
    
    data_bits = generate_random_bits(packet_size)
    actual_channel = channel_vals[current_index:current_index + packet_size]
    received_data = transmit_data(data_bits, actual_channel, snr_db)
    
    dbar = np.sign(np.real(received_data * np.conjugate(htilde) / (np.abs(htilde)**2)))
    hnought = estimate_channel(received_data, dbar)
    mse = np.abs(htilde - hnought)**2
    
    if mse > mse_threshold:
        new_estimates, current_index = collect_pilot_estimates(current_index, initial_pilots)
        channel_estimates.extend(new_estimates)
        pilot_requests.extend(range(current_index - len(new_estimates) * pilot_size, current_index, pilot_size))
        total_bits += len(new_estimates) * pilot_size
    else:
        channel_estimates.append(hnought)
        current_index += packet_size
        total_bits += packet_size

pilot_ratio = 1/(total_bits / (len(pilot_requests) * pilot_size))
print(f"{pilot_ratio}")
END
)
        
        # Append results to output file
        echo "$dopp,$snr,$ratio" >> $output_file
        
    done
done
# Create the plot using Python
python3 - <<END
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Read the data
df = pd.read_csv("$output_file")

# Create the plot
plt.figure(figsize=(12, 8))

# Plot a line for each SNR value
for snr in df['SNR'].unique():
    data = df[df['SNR'] == snr]
    plt.plot(data['dopp'], data['ratio'], label=f'SNR={snr:.3f}', marker='o', markersize=4)

plt.title('Pilot Ratio vs SNR')
plt.xlabel('Doppler ratio')
plt.ylabel('Pilot Ratio')
plt.grid(True)
plt.legend(loc='upper right')
plt.tight_layout()

# Save the plot
plt.savefig("ratio_vs_SNR_plot.png", bbox_inches='tight', dpi=300)
plt.close()

print("Created plot: ratio_vs_SNR_plot.png")
END