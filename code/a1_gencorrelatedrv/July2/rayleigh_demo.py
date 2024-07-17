import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, rayleigh

N = 100000

def generate(n_samples=N):
    gaussian1 = np.random.normal(0, 1, n_samples)
    gaussian2 = np.random.normal(0, 1, n_samples)

    rayleigh_rv = np.sqrt(gaussian1**2 + gaussian2**2)

    # Plot PDFs
    fig, ax = plt.subplots(2, 2, figsize=(14, 12))

    # Plot Gaussian PDFs
    x = np.linspace(-4, 4, 1000)
    ax[0, 0].plot(x, norm.pdf(x, 0, 1), label='Gaussian PDF', color='blue')
    ax[0, 0].hist(gaussian1, bins=100, density=True, alpha=0.6, color='blue', edgecolor='black')
    ax[0, 0].hist(gaussian2, bins=100, density=True, alpha=0.6, color='red', edgecolor='black')
    ax[0, 0].set_title('Gaussian Random Variables')
    ax[0, 0].set_xlabel('Value')
    ax[0, 0].set_ylabel('Density')
    ax[0, 0].legend()

    # Plot Rayleigh PDF
    x = np.linspace(0, 6, 1000)
    ax[0, 1].plot(x, rayleigh.pdf(x, scale=1), label='Rayleigh PDF', color='green')
    ax[0, 1].hist(rayleigh_rv, bins=100, density=True, alpha=0.6, color='green', edgecolor='black')
    ax[0, 1].set_title('Rayleigh Random Variable')
    ax[0, 1].set_xlabel('Value')
    ax[0, 1].set_ylabel('Density')
    ax[0, 1].legend()

    # Plot Gaussian random variables as series
    ax[1, 0].plot(gaussian1[:1000], label='Gaussian RV 1', color='blue')
    ax[1, 0].plot(gaussian2[:1000], label='Gaussian RV 2', color='red')
    ax[1, 0].set_title('Gaussian Random Variables Series')
    ax[1, 0].set_xlabel('Sample Index')
    ax[1, 0].set_ylabel('Value')
    ax[1, 0].legend()

    # Plot Rayleigh random variable as series
    ax[1, 1].plot(rayleigh_rv[:1000], label='Rayleigh RV', color='green')
    ax[1, 1].set_title('Rayleigh Random Variable Series')
    ax[1, 1].set_xlabel('Sample Index')
    ax[1, 1].set_ylabel('Value')
    ax[1, 1].legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    generate()

