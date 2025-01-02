import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
import json

# Function to generate Rayleigh fading data
def generate_rayleigh_fading_data(num_samples):
    complex_data = (1/np.sqrt(2)) * (np.random.randn(N) + 1j * np.random.randn(N))
    return complex_data

# Define parameters
N = 500000  # Number of samples

#Simulation of BPSK BER with generated channel data
with open("channels.txt", "r") as f:
    channel_vals = [json.loads(line)[0] for line in f]
    channel_vals = np.array([complex(re, im) for re, im in channel_vals])
print(channel_vals[:20])
ip = np.random.rand(N) > 0.5
s = 2 * ip - 1

Eb_N0_dB = np.arange(-3, 36)
nErr = np.zeros(len(Eb_N0_dB))

for i in range(len(Eb_N0_dB)):
    noise = (1/np.sqrt(2)) * (np.random.randn(N) + 1j * np.random.randn(N))
    h = channel_vals  # Use the loaded channel values
    y = h * s + 10**(-Eb_N0_dB[i] / 20) * noise
    yHat = y * np.conj(h) / np.abs(h) ** 2
    ipHat = np.real(yHat) > 0
    nErr[i] = np.sum(ip != ipHat)

simBer = nErr / N
theoryBerAWGN = 0.5 * erfc(np.sqrt(10**(Eb_N0_dB / 10)))
EbN0Lin = 10**(Eb_N0_dB / 10)
theoryBerRayleigh = 0.5 * (1 - np.sqrt(EbN0Lin / (EbN0Lin + 1)))

plt.figure()
plt.semilogy(Eb_N0_dB, theoryBerAWGN, 'cd-', linewidth=2, label='AWGN Theory')
plt.semilogy(Eb_N0_dB, theoryBerRayleigh, 'bp-', linewidth=2, label='Rayleigh Theory')
plt.semilogy(Eb_N0_dB, simBer, 'mx-', linewidth=2, label='Simulation')

plt.axis([-3, 35, 1e-5, 0.5])
plt.grid(True, which='both')
plt.legend()
plt.xlabel('Eb/No (dB)')
plt.ylabel('Bit Error Rate (BER)')
plt.title('BER for BPSK in Rayleigh Fading Channel')
plt.show()
