import numpy as np
import matplotlib.pyplot as plt
from statsmodels.tsa.arima_process import ArmaProcess

# Given autocorrelation function R(k)
R = [1, 0.6, 0.36, 0.216]  # Example values, this should be your actual R(k)

# Solve Yule-Walker equations to find AR coefficients
def yule_walker(R, order):
    from scipy.linalg import toeplitz
    r = np.array(R[:order])
    R_matrix = toeplitz(r)
    phi = np.linalg.inv(R_matrix).dot(r)
    return phi

# Estimate AR(1) coefficient
phi = yule_walker(R, 1)[0]

# Assume an MA(1) process, estimate theta based on R(1)
theta = 0.5  # This is a placeholder; more complex methods are required for exact estimation

# ARMA(1,1) coefficients
ar_params = [1, -phi]
ma_params = [1, theta]

# Generate ARMA(1,1) process
arma_process = ArmaProcess(ar_params, ma_params)
np.random.seed(42)
X = arma_process.generate_sample(nsample=1000)

# Plot the time series
plt.figure(figsize=(10, 6))
plt.plot(X, marker='o', markersize=2)
plt.title('Generated ARMA(1,1) Process')
plt.xlabel('Time')
plt.ylabel('Value')
plt.grid(True)
plt.tight_layout()
plt.show()

# Plot autocorrelation
from statsmodels.graphics.tsaplots import plot_acf
plot_acf(X, lags=20)
plt.show()
