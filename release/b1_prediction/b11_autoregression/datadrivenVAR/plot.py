import numpy as np
import matplotlib.pyplot as plt

results = np.loadtxt("results.txt", delimiter=",", skiprows=1)
pilot_ratios = results[:, 3]
mses = results[:, 4]
bers = results[:, 5]

plt.figure(figsize=(12, 8))
plt.subplot(3, 1, 1)
plt.plot(pilot_ratios, label="Pilot Ratios")
plt.ylabel("Pilot Ratio")
plt.legend()

plt.subplot(3, 1, 2)
plt.plot(mses, label="MSE")
plt.ylabel("MSE")
plt.legend()

plt.subplot(3, 1, 3)
plt.plot(bers, label="BER")
plt.ylabel("BER")
plt.xlabel("Configurations")
plt.legend()

plt.tight_layout()
plt.show()
