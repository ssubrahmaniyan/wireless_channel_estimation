# Channel Generation and Estimation
### Concepts covered
Any wireless communication medium is characterized by a channel gain, which is a complex value that is present as a multiplicative factor on the symbols being transmitted through the channel. If we wish to 
transmit a bit of information, with the BPSK modulation scheme, the symbol x that we transmit will be either 1 (in case the bit we want to communicate is 0) or -1 (in case the bit we want to communicate is 1). 
Then, the symbol y that is received is modeled as: 

<div align="center">
    y = h*x + n
</div>

where n is noise of a specified distribution and level determined by the Signal-to-Noise Ratio (SNR) of the transmission. 

The symbol y is converted back into a received bit, based on the modulation scheme we choose to use. In case of BPSK, the sign of the estimate of x, given a value for y, is used to conclude the nature of the bit
that was trasmitted. If the estimate of x is positive, then the bit that was transmitted was 0, and if the estimate of x is negative, then the bit that was transmitted was 1. 

A crucial step in estimating x given a value for y is the determination of the channel value to be used. One approach, which has been used in the plots below, is to assume we know the actual value of the channel for every bit being transmitted. Another approach involves transmitting headers, or 
packets of information where we know the exact bits in the packets, and then determining a value of the channel for these header bits. We then assume the channel values so determined are constant over the next few bits we transmit, and use these channel values to demodulate and understnad the next few data bits.
Other approaches that reduce the involvement of headers are being explored as a part of the project. 

A metric for evaluating whether our estimate of the bits transmitted is accurate is compute the BER of the predicted signal against the actual signal transmitted. Since this value would nturally improve as the SNR improves, a plot of the BER Vs SNR for the channel model and the modulation scheme would be informative. 
This plot has been made below. 

3 types of channel values I have examined in the plot below are:
* AWGN channels
* Uncorrelated Rayleigh faded channels with AWGN noise
* Correlated Rayleigh faded channels with AWGN noise

Channel impairments that have been covered as a part of the channel models explored above are:
* Fading: Variations in signal amplitude due to multipath.
* Doppler Shift: Frequency change due to relative motion between transmitter and receiver.
* Noise: Random additive disturbance modeled as Gaussian noise.

Correlation in the Rayleigh Fading channel has been introduced by the sum of sinusoids method as outlined in <https://www.ee.iitm.ac.in/~giri/pdfs/EE5141/Jakes-Simulation.pdf>. This has been adapted from the Rapport textbook and can be easily implemented to generate a Jakes channel.

### Plots
![plot](/Release/Channel Generation and Estimation/BER_Vs_SNR.png)

This plot shows the BER vs SNR characteristics for the 3 channels described above. Clearly, an AWGN channel is the most SNR dependent, and in the presence of SNR values of around 10 dB, the transmitted bits will be perfectly received. 
Rayleigh fading poses a slightly tougher problem of bit estimation. 

Need to still understand why the Jakes spectrum plot is coinciding with the uncorrelated Rayleigh faded channel plot. 


