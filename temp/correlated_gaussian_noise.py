import numpy as np
import matplotlib.pyplot as plt

# Given parameters
K1 = 1.5
K2 = 2.0
K_12 = 1.2
rho = K_12 / np.sqrt(K1 * K2)
Delta_t = 0.01
n = 1000

# Generate random normal variables
xi1 = np.random.normal(0, 1, n)
xi2 = np.random.normal(0, 1, n)

# Generate the two correlated processes
W1 = np.sqrt(2 * np.pi * K1 / Delta_t) * xi1
W2_prime = np.sqrt(2 * np.pi * K2 / Delta_t) * xi2
W2 = rho * W1 + np.sqrt(1 - rho**2) * W2_prime

# Time vector
t = np.arange(n) * Delta_t

# Plot the two processes
plt.figure(figsize=(12, 6))
plt.plot(t, W1, label='W1')
plt.plot(t, W2, alpha=0.75, label='W2')
plt.title('Two Correlated Processes')
plt.xlabel('Time')
plt.ylabel('Value')
plt.legend()
plt.show()

# Compute rolling correlation
window_size = 500
correlation = [np.corrcoef(W1[i:i+window_size], W2[i:i+window_size])[0, 1] for i in range(n - window_size + 1)]
t_corr = t[window_size-1:]

# Plot the correlation
plt.figure(figsize=(12, 6))
plt.plot(t_corr, correlation, label='Rolling Correlation')
plt.title('Rolling Correlation between W1 and W2')
plt.xlabel('Time')
plt.ylabel('Correlation')
plt.ylim(-1, 1)
plt.legend()
plt.show()


