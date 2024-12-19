# Data transmissions modelling

## Theoretical estimates of channel models

The theoreticalplots.py file can be used to make the theoretical plots for BPSK with different channels. It shows the expected behaviour using the determining equations and results of a simulation in a Rayleigh based channel.

![plot](/release/a2_transmission_modelling/fig0_theoreticalplots.png)

## Modelling when all the channels are known

The file allchannelknown.py models transmission in a fading channel with the assumption that the channel values are known at each instant and computes the BER with differing values of the SNR of the signal. 

The file should be run by generating channel values of the required doppler frequency and modifying the doppler value in the generating function at the end of the allchannelknown.py file. 

    return acorr / zero_lag

Two representative plots are as follows:

![plot](/release/a2_transmission_modelling/fig1_allchannelknown10hz.png)

This plot shows the comparison between the theoretical values of the BER as expected from a simple Rayleigh fading channel and the BER for a channel modelled with the Jakes spectrum with a doppler shift of 10 Hz. As is expected, the error rates are higher in the doppler shifted channel due to the variation in the channel itself as a function of time

![plot](/release/a2_transmission_modelling/fig2_allchannelknown100hz.png)

Similarly, this plot shows the BER for a channel with a doppler shift of 100Hz. **Is it sensible that they look very similar?**

## Full transmission of data
        alpha_m = np.array([((2 * np.pi * n) - np.pi + al)/(4*N) for n, al in enumerate(alpha)])
