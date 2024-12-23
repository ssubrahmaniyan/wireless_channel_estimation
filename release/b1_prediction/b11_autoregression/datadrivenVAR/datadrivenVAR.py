import numpy as np
import matplotlib.pyplot as plt
from statsmodels.tsa.api import VAR
import json

with open("channels.txt", "r") as f:
    channel_vals = [json.loads(line)[0] for line in f]
channel_vals = np.array([complex(re, im) for re, im in channel_vals])

Fd = 100
Fs = 100000
quanta = 5
pilot_size = 50
data_size = quanta * pilot_size
packet_size = data_size + pilot_size
var_order = 25
mse_threshold = 0.1
snr_db = 10
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

def find_ber(decoded_data, original_data):
    count = 0
    for i in range(len(decoded_data)):
        if decoded_data[i] != original_data[i]:
            count += 1
    return count/len(decoded_data)


current_index = 0
channel_estimates = []
pilot_requests = []
total_bits = 0

initial_estimates, current_index = collect_pilot_estimates(current_index, initial_pilots)
channel_estimates.extend(initial_estimates)
pilot_requests.extend(range(0, current_index, pilot_size))
total_bits += initial_pilots * pilot_size

headermses = []
databers = []
if len(channel_estimates) <= var_order:
    raise ValueError("Insufficient data for initial estimations")

while current_index < len(channel_vals) - packet_size:
    try:
        htilde = predict_channel_var(channel_estimates)
    except ValueError as e:
        break
    
    pilot_bits = generate_random_bits(pilot_size)
    data_bits = generate_random_bits(data_size)
    packet = np.concatenate((pilot_bits, data_bits))
    actual_channel = channel_vals[current_index:current_index + packet_size]
    received_data = transmit_data(packet, actual_channel, snr_db)
    pilot_requests.append(1)
    # plt.plot(packet)
    # plt.plot(received_data.real)
    # plt.plot(received_data.imag)
    # plt.show()
    # Estimate the channel from the header bits and add it to the VAR model
    hbar = estimate_channel(received_data[:pilot_size], pilot_bits)
    # Now the first part of the estimation is done
    # Demodulate a part of the data bits
    channel_estimates.append(hbar)
    hestimates = [hbar]
    hpred = [htilde]
    bers = []
    for i in range(quanta):
        try:
            htilde = predict_channel_var(channel_estimates)
        except ValueError as e:
            break
        
        hpred.append(htilde)
        dbar = np.sign(np.real(received_data[(1 + i)*pilot_size:(2 + i) * pilot_size] * np.conjugate(htilde) / (np.abs(htilde)**2)))
        bers.append(find_ber(dbar, packet[(1 + i)*pilot_size:(2 + i) * pilot_size]))
        hnought = estimate_channel(received_data[(1 + i)*pilot_size:(2 + i) * pilot_size], dbar)
        channel_estimates.append(hnought)
        hestimates.append(hnought)
        # print(dbar)
        # print(packet[(1 + i)*pilot_size:(2 + i) * pilot_size])
        # print(hnought, hbar)
    

    total_bits += packet_size
    current_index += packet_size
    hpred = np.array(hpred)
    hestimates = np.array(hestimates)
    mse = np.sqrt(np.sum(np.abs(hpred - hestimates)**2)/len(hpred))
    databers.append(np.sum(bers)/quanta)
    headermses.append(mse)

    if mse > mse_threshold:
        new_estimates, current_index = collect_pilot_estimates(current_index, initial_pilots)
        channel_estimates.extend(new_estimates)
        pilot_requests.extend(range(current_index - len(new_estimates) * pilot_size, current_index, pilot_size))
        total_bits += initial_pilots * pilot_size
    
pilot_ratio = (len(pilot_requests) * pilot_size/total_bits)
print(f"Pilot ratio = {pilot_ratio}")

print(f"Average MSE = {np.mean(headermses)}")
print(f"Average BERS = {np.mean(databers)}")