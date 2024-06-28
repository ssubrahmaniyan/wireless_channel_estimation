import numpy as np
import matplotlib.pyplot as plt
from statsmodels.tsa.arima_process import ArmaProcess

phi = 0.6
theta = 0.4
n = 1000

ar = np.array([1, phi])
ma = np.array([1, theta])
arma_process = ArmaProcess(ar, ma)

uncorrelated_variable = np.random.normal(0, 1, n)
correlated_variable = arma_process.generate_sample(n) + uncorrelated_variable

plt.figure(figsize=(12, 6))
plt.plot(uncorrelated_variable, color = 'r')
plt.plot(correlated_variable, alpha=0.75, color = 'b')
plt.show()

