import numpy as np
import matplotlib.pyplot as plt
from statsmodels.tsa.api import VAR
import json

# Load the complex channel values
with open("channels.txt", "r") as f:
    channel_vals = [json.loads(line)[0] for line in f]
channel_vals = np.array([complex(re, im) for re, im in channel_vals])

# Parameters
packet_size = 50
pilot_size = 50
var_order = 25
mse_threshold = 0.1
snr_db = 10  # SNR in dB
initial_pilots = var_order + 2  # Minimum number of samples needed for VAR fitting + buffer

def generate_random_bits(n):
    return np.random.choice([-1, 1], size=n)

def transmit_data(data, channel, snr_db):
    signal = channel * data
    signal_power = np.mean(np.abs(signal)**2)
    noise_power = signal_power / (10**(snr_db/10))
    noise = np.sqrt(noise_power/2) * (np.random.randn(len(data)) + 1j * np.random.randn(len(data)))
    return signal + noise

def estimate_channel(received, data):
    """
    Estimate a single constant channel value for the entire packet using least squares
    """
    X = data.reshape(-1, 1)
    h = np.linalg.lstsq(X, received, rcond=None)[0][0]
    return h

def collect_pilot_estimates(current_index, num_pilots):
    """
    Collect channel estimates from multiple pilot transmissions
    """
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
    """
    Train VAR model on channel history and predict next channel value
    """
    if len(channel_history) <= var_order:
        raise ValueError(f"Not enough samples for VAR. Need >{var_order}, got {len(channel_history)}")
        
    channel_history_separated = np.column_stack((np.real(channel_history), np.imag(channel_history)))
    model = VAR(channel_history_separated)
    results = model.fit(var_order)
    forecast = results.forecast(channel_history_separated[-var_order:], steps=steps)
    return forecast[0, 0] + 1j * forecast[0, 1]  # Return only next value

# Initialize
current_index = 0
channel_estimates = []  # Store all channel estimates (both from pilots and data)
pilot_requests = []
total_bits = 0

# Initial pilot transmissions - ensure enough samples for VAR
initial_estimates, current_index = collect_pilot_estimates(current_index, initial_pilots)
channel_estimates.extend(initial_estimates)
pilot_requests.extend(range(0, current_index, pilot_size))
total_bits += initial_pilots * pilot_size

if len(channel_estimates) <= var_order:
    raise ValueError(f"Could not collect enough initial pilots. Need >{var_order}, got {len(channel_estimates)}")

# Main loop
while current_index < len(channel_vals) - packet_size:
    # Predict channel using VAR
    try:
        htilde = predict_channel_var(channel_estimates)
    except ValueError as e:
        break
    
    # Data transmission
    data_bits = generate_random_bits(packet_size)
    actual_channel = channel_vals[current_index:current_index + packet_size]
    received_data = transmit_data(data_bits, actual_channel, snr_db)
    
    # Demodulate using predicted channel
    dbar = np.sign(np.real(received_data * np.conjugate(htilde) / (np.abs(htilde)**2)))
    
    # Re-estimate channel using demodulated data
    hnought = estimate_channel(received_data, dbar)
    
    # Calculate error between predicted and re-estimated channel
    mse = np.abs(htilde - hnought)**2
    
    if mse > mse_threshold:
        # Send new pilots if error is too large
        new_estimates, current_index = collect_pilot_estimates(current_index, initial_pilots)
        channel_estimates.extend(new_estimates)
        pilot_requests.extend(range(current_index - len(new_estimates) * pilot_size, current_index, pilot_size))
        total_bits += len(new_estimates) * pilot_size
    else:
        # Use hnought as new channel estimate and continue
        channel_estimates.append(hnought)
        current_index += packet_size
        total_bits += packet_size

# Calculate pilot efficiency
pilot_ratio = 1/(total_bits / (len(pilot_requests) * pilot_size))
print(pilot_ratio)