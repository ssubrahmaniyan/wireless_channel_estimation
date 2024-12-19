import json
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import jn  # Bessel function of the first kind

def autocorrelation(x):
    result = np.correlate(x, x.conjugate(), mode='full')
    return result

def normalize_autocorrelation(acorr):
    # Normalize by the zero lag value (maximum value)
    zero_lag = acorr[len(acorr) // 2]
    return acorr / zero_lag

# Load and process the channel data
with open("channels.txt", "r") as f:
    # Deserialize each line into complex numbers
    channel_values = [
        np.array([complex(re, im) for re, im in json.loads(line)])
        for line in f
    ]

# Convert the list of NumPy arrays to a 2D NumPy array
channel_values = np.array(channel_values)
channel_count = channel_values.shape[1]

# Time spacing between samples
time_spacing = 1 / 100000

# Prepare to plot
plt.figure(figsize=(14, 6 * channel_count))  # Adjust figure size to fit both plots

for i in range(channel_count):
    # Extract channel data for plotting
    channel_data = channel_values[:100000, i]
    
    # Compute the logarithmic amplitude for the signal
    signal_amplitude = 10 * np.log10(np.abs(channel_data))
    
    # Compute the autocorrelation of the signal
    channel_autocorr = autocorrelation(channel_data)
    
    # Normalize the autocorrelation
    normalized_autocorr = normalize_autocorrelation(channel_autocorr)
    
    # Center around zero
    lags = np.arange(-len(channel_data) + 1, len(channel_data))
    mid = len(lags) // 2
    centered_autocorr = normalized_autocorr[mid - len(lags)//2 : mid + len(lags)//2 + 1]
    centered_lags = lags[mid - len(lags)//2 : mid + len(lags)//2 + 1]

    # Compute J0(200 * pi * t) where t is the lag in seconds
    bessel_func = jn(0, 200 * np.pi * centered_lags * time_spacing)

    # Create subplots: one for the signal and one for the autocorrelation
    plt.subplot(channel_count, 2, i * 2 + 1)
    plt.plot(np.arange(len(signal_amplitude)), signal_amplitude)
    plt.title(f'Channel {i+1} - Signal Amplitude (dB)')
    plt.xlabel('Sample Index')
    plt.ylabel('Amplitude (dB)')

    plt.subplot(channel_count, 2, i * 2 + 2)
    plt.plot(centered_lags * time_spacing, (centered_autocorr), label='Normalized Autocorrelation')
    plt.plot(centered_lags * time_spacing, bessel_func, label='J0(200πt)', linestyle='--')
    plt.title(f'Channel {i+1} - Normalized Autocorrelation and Bessel Function')
    plt.xlim(-0.05, 0.05)
    plt.xlabel('Lag (seconds)')
    plt.ylabel('Value')
    plt.legend()

    # Find the index where the autocorrelation first reaches 1
    first_one_index = np.where(centered_autocorr == 1)[0]
    
    if first_one_index.size > 0:
        first_one_index = first_one_index[0]  # Get the first index
        print(f'Channel {i+1}: Autocorrelation first reaches 1 at lag index {first_one_index} (lag: {centered_lags[first_one_index] * time_spacing:.6f} seconds)')

        # Find the first index after first_one_index where the autocorrelation drops below 0.7
        subsequent_drop_index = np.where(centered_autocorr[first_one_index:] < 0.7)[0]
        if subsequent_drop_index.size > 0:
            drop_index = subsequent_drop_index[0] + first_one_index  # Adjust index based on the slicing
            print(f'Channel {i+1}: Autocorrelation drops below 0.7 at lag index {drop_index} (lag: {centered_lags[drop_index] * time_spacing:.6f} seconds)')
        else:
            print(f'Channel {i+1}: Autocorrelation does not drop below 0.7 after reaching 1.')

# Show the plots
plt.tight_layout()
plt.show()

