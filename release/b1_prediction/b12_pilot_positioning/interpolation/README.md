# Doppler filter based channel generation

Firstly, to fix the rogue channel generation program, I plotted the filter response to see if I was making a mistake there.
![plot](/filter_spectrum.png)

Because the spectrum looks as it should be, I went on to checking a few other metrics to see how they correlated with the outputs that I was seeing. To do this, I made a plot of the distribution of the magnitudes and phases of the channel points that I had generated using the filter. I was able to infer that a distribution that is similar to the Rayleigh's distribution would give an BER curve similar to the Rayleigh and so on. Two relevant samples are below.

### 100 Hz
![plot](/distro_100hz.png)

The corresponding transmission plot(when each point's channel is known):

![plot](/transmit_100hz.png)

For the same number of data points and sampling frequency, but with a lower doppler,

### 10 Hz

![plot](/distro_10hz.png)

The corresponding transmission plot is:

![plot](/transmit_10hz.png)

### 50 Hz

![plot](/distro_50hz.png)

The corresponding transmission plot is:

![plot](/transmit_50hz.png)


## Observations

1. I noticed that the plot very quickly tends towards either of the edge cases(AWGN and Rayleigh) instead of a slow progression. This is very clear in the above case where while 1Hz is practically the same as AWGN, 100Hz is the same as Rayleigh. 
1. From the previous observation, I attempted to vary the Fs for a fixed Fd. The observation was as expected, with an increase in Fs quickly bringing the curve downwards towards the AWGN curve and a decrease doing the opposite. This means Fd/Fs is more relevant than Fd itself.
1. The distribution curves given above do not look smooth at all. I was able to make them smoother by generating and averaging out more. This also solved my problem of slightly inconsistent results(where the transmission plots for the same parameters would look rather different). I think this means 10^6 points is not enough for the models and 10^8 might be a better option. But is the non-smoothness important as being smooth would mean a behavior similar to that of the Rayleigh spectrum? A sampling frequency of about 10^8 also looks like a good option.

I have used the filter method for channel generation and not the sum of sinusoids for the following results.

# Simple nointerpolation comparison for different doppler

A simple transmission method where the channel is estimated from the pilots and then assumed constant for the data part of the packet. The assumed channel is used to demodulate the data as expected. 

