# Error Correcting Codes
### Algorithm

The process of receiving symbols after being transmitted through the channel and then demodulating them to predict a set of bits needs to be accurate if we are to get a reliable estimate of the channel. The demodulation process may induce some errors, and a more reliable method of ensuring that the demodulation being performed is accurate is if we use Error Correcting Codes that introduce parity bits into the message we are transmitting. 
A popular method of ECC for wireless communication in the 5G standard is LDPC - Low Density Parity Check Codes. I have observed implementations of encoding and decoding messages using LDPC with a library named pyldpc. 
A few shortcomings of the implementation in the library have also been fixed and integrated in my codes, which can be found in LDPC.ipynb, along with all the code used for plotting. 

The general algorithm involving LDPC is given:
* Generate a set of bits to be transmitted
* Encode the bits into LDPC format and convert the bits to symbols as per the modulation scheme (BPSK here)
* Transmit the bits through the channel and receive them
* Demodulate the entire received block of symbols into bits
* Decode the LDPC code back into the recevied message and compute BER on the message

LDPC Encoding -
* using a sparse parity-check matrix to place constraints on the parity bits
* uses a generator matrix to actually get the parity bits from the data bits
LDPC Decoding - for our channel, uses the sum-product message passing/belief propagation algorithm - need to understand in more detail

Parameters:
- k = the number of bits to be transmitted originally before coding - this is the data we want to send per packet, so more packet size increases latency
- n = the number of bits to be transmitted after coding is done on k bits, so n-k parity bits is used
- H = parity check matrix
- G = generator matrix
- d_v = degree of the variable node - number of 1s per column in H
- d_c = degree of the check node - number of 1s per row in H
- d_v and d_C along with n specify the structure of H and G.
- H has n columns and (n * d_v/d_c) = m rows.
- The rate of coding R = 1 - (m/n) = k/n
- Closer R is to 1, more uncoded the data is, suitable for less noise or no noise.

The effect of coding on different channels, provided we know the channel, has been examined below. One issue is that in the current implementation, k/n is not exaclty equal to the rate. This issue needs to be fixed. 

### Plots
![BER Vs SNR Plot1](/Release/Error%20Correction%20Codes/BER_Vs_Eb_AWGN_Coding.png)
The plot above shows the BER Vs SNR characteristic for an AWGN channel where we assume we know the channel.  

![BER Vs SNR Plot2](/Release/Error%20Correction%20Codes/BER_Vs_Eb_Ray_Coding.png)
The plot above shows the BER Vs SNR characteristic for an uncorrelated Rayleigh channel where we assume we know the channel.

![BER Vs SNR Plot3](/Release/Error%20Correction%20Codes/BER_Vs_Eb_Jakes_Coding.png)
The plot above shows the BER Vs SNR characteristic for a correlated Rayleigh channel where we assume we know the channel.
