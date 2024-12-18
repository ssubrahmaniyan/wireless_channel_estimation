# Wireless_Channel_Estimation

Workflow for SISO wireless communication:
- Generate a signal to be transmitted, with some bits in each packet as reference bits
- Receive the modified packet at the receiver
- Denoise the packet using least squares estimation, y = Ax, find an optimal x so MSE(y, Ax) is minimised
- Use the reference bits and determine the channel for the packet (assumed to be time-invariant for the duration of the packet)
- Use this channel to then infer the remaining of the original signal from the modified packet, iterate over whole signal

# Concepts examined as of 15/12/2024

- Have understood how to generate 1D real or complex channels that have the same autocorrelation as Jakes' Sum of Sinusoids method with Rayleigh Multipath Fading
- Have implemented simple binary signal generation, transmission and BPSK demodulation
- Have implemented the standard header-packet channel estimation pipeline, and observed BER vs SNR (waterfall) plots and varied header:packet ratios too
- Have implemented a Vector Autoregression based model to predict channel gains into the future after training on least-squares based channels developed by sending pilots where data is known.
- Have implemented 2 autoregressors - one for magnitude and phase each - instead of a VAR model, and observed that the prediction accuracy is not as good since VAR maintains the autocorrelation better.
- Plotted retransmission frequency against doppler frequency for same and varying MSE thresholds. 
