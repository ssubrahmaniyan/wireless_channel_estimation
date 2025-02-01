# Error Correcting Codes
### Algorithm

The process of receiving symbols after being transmitted through the channel and then demodulating them to predict a set of bits needs to be accurate if we are to get a reliable estimate of the channel. The demodulation process may induce some errors, and a more reliable method of ensuring that the demodulation being performed is accurate is if we use Error Correcting Codes that introduce parity bits into the message we are transmitting. 
A popular method of ECC for wireless communication in the 5G standard is LDPC - Low Density Parity Check Codes. 

I have observed implementations of encoding and decoding messages using LDPC with a library named pyldpc.  A few shortcomings of the implementation in the library have also been fixed and integrated in my codes, which can be found in LDPC.ipynb, along with all the code used for plotting. 

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

The effect of coding on different channels, provided we know the channel, has been examined below. Note that the rate is defined by k/n, and not by d_v/d_c. 

### Plots
<img src="/Release/Error%20Correction%20Codes/BER_Vs_Eb_AWGN_Coding.png" width="49%">

The plot above shows the BER Vs SNR characteristic for an AWGN channel where we assume we know the channel.  


<img src="/Release/Error%20Correction%20Codes/BER_Vs_Eb_Ray_Coding.png" width="49%">

The plot above shows the BER Vs SNR characteristic for an uncorrelated Rayleigh channel where we assume we know the channel.


<img src="/Release/Error%20Correction%20Codes/BER_Vs_Eb_Jakes_Coding.png" width="49%">

The plot above shows the BER Vs SNR characteristic for a correlated Rayleigh channel where we assume we know the channel, for 100000 bits and Fs = 100000.


<img src="/Release/Error%20Correction%20Codes/BER_Vs_Eb_Jakes_Coding2.png" width="49%">

The plot above shows the BER Vs SNR characteristic for a correlated Rayleigh channel where we assume we know the channel, for 10000 bits and Fs = 10000.


<img src="/Release/Error%20Correction%20Codes/BER_Vs_SNR.png" width="49%">

For reference, the plot above shows the BER Vs SNR characteristic for 3 different channel types without coding. The results suggest that the Jakes method being suggested is very similar to a typical Rayleigh channel generator. This needs to be discussed. 


### Encoding Window Length for correlated channels
The encoding window length (n) that is used in case of correlated channels must be chosen carefully with respect to the coherence time - the time for which the channel can be assumed to be constant. This time depends inversely with the Doppler Frequency and this relation can be assumed to be an exact inverse. 

Hence for a Doppler frequency of 100 Hz, the coherence time is 0.01 seconds, and for a Sampling frequency of 10000 Hz, this means that 100 datapoints are transmitted each second. Ideally we would want 1 in 4 coherence time packets to be entirely corrected by LDPC, so in 400 samples that are transmitted, we can correct upto 100 errors. 

In order to correct t errors, we need the minimum distance of the code to satisfy this:
<div align="center">
    d_min > 2*t + 1
</div>

The minimum distance is proportional to the number of parity bits, and hence as we increase n, we can correct more errors in the transmission and achieve lower BER at the same SNR. For a rate of 0.5, we need to transmit 400 bits or 800 symbols with 400 parity bits to achieve a redundancy of 100 bits in the original message. The simulation below explores the BER at different n values or encoding lengths. 

The plots below show the variation in coding gain for different encoding window lengths for 10000 bits transmitted. 

### Plots for variable encoding window length

<img src="/Release/Error%20Correction%20Codes/SNR_Vs_Eb_AWGN_Nvals.png" width="49%">

The plot above shows the BER vs SNR plots for different encoding window lengths for 10000 bits and Fs = 10000 for a AWGN channel

<img src="/Release/Error%20Correction%20Codes/SNR_Vs_Eb_Jakes_Nvals.png" width="49%">

The plot above shows the BER vs SNR plots for different encoding window lengths for 10000 bits and Fs = 10000 for a correlated Jakes channel

### Hard and Soft Decision Decoding and accounting for the Fd/Fs ratio
<img src="/Release/Error%20Correction%20Codes/BER_Vs_EB_3_channel.png" width="49%">

The plot above shows the BER vs SNR plots for different channels for Fd = 1 Hz and Fs = 1000000 Hz for 10000 points. Clearly, the Jakes curve lies between the AWGN and the Rayleigh curves at appropriate values of Fd. 

The algorithm used so far has involved digitizing the received signals, after dividing the received symbols with the known or estimated channel, and predicting the transmitted bits from this digitization. Such an approach is called a Hard Decision Decoder. 

To improve the accuracy of the prediction, a Soft Decision Decoder is used. There ae different algorithms through which this may be done as specified by https://ieeexplore.ieee.org/document/7050112 and https://iopscience.iop.org/article/10.1088/1757-899X/1105/1/012039.

Approaches tried out:
- Computing Log Likelihood Ratio (LLR) for received symbols after normalizing with the channel. For BPSK, the LLR is given by
 <div align="center">
    LLR = 2 * y / variance of noise
</div> 
LLR > 0 implies that the bit is 0, and LLR < 0 implies that the bit is 1. This approach is practically the same as a hard decision decoder. 

<img src="/Release/Error%20Correction%20Codes/BER_Vs_Eb_Ray_LLR.png" width="49%">

The plot above shows the BER vs Eb/N0 plot for Rayleigh fading channel for 10000 bits with this approach. 


- Computing sum of squared error between normalized recieved symbols and all possible transmitted symbols whose encoded versions have the least MSE with the received normalized symbol vector - highly inefficient and does not scale with k. The paper cited above only suggests a way to speed it up on hardware. 

<img src="/Release/Error%20Correction%20Codes/BER_Vs_Eb_Ray_MSE.png" width="49%">

The plot above shows the BER vs Eb/N0 plot for Rayleigh fading channel for 10000 bits with this approach. 


- LDPC decoding with LLRs - the normalized received symbols are used in the Belief Propagation Algorithm and the posterior LLR signs are used to predict the transmitted bits.


<img src="/Release/Error%20Correction%20Codes/BER_Vs_Eb_AWGN_LogBP.png" width="49%">

The plot above shows the BER vs Eb/N0 plot for AWGN channel with this approach

# LDPC with Encoding and Soft-Input Iterative Decoding using LLRs

<img src="/Release/Error%20Correction%20Codes/AWGN_Final2.png" width="48%"> <img src="/Release/Error%20Correction%20Codes/Reference_AWGN.png" width="42%">

The plot above on the left shows the BER vs Eb/N0 plot for an AWGN channel with a rate of 0.5 and codeword length of 648. The plot above on the right shows the expected BER vs Eb/N0 behaviour specified in the references - Error Correction Coding by Todd K. Moon. 


<img src="/Release/Error%20Correction%20Codes/AWGN_Multiple_Lengths.png" width="60%">

The plot above shows the BER vs Eb/N0 plot for different rates and codeword lengths used for an AWGN channel. 


<img src="/Release/Error%20Correction%20Codes/Rayleigh_Final.png" width="49%"> <img src="/Release/Error%20Correction%20Codes/Reference_Rayleigh.png" width="41%">

The plot above on the left shows the BER vs Eb/N0 plot for an uncorrelated Rayleigh channel with a rate of 0.5 and codeword length of 648. The plot on the right shows the expected Rayleigh plot from the references. 


<img src="/Release/Error%20Correction%20Codes/Rayleigh_Ideal.png" width="49%">

The plot above shows the BER vs Eb/N0 plot for an uncorrelated Rayleigh channel with rate of 0.50 and codeword length of 2040

<img src="/Release/Error%20Correction%20Codes/BER_Vs_EB_3_channel_50Hz.png" width="49%">

For reference, the plot above shows the BER vs Eb/N0 plot for 3 different channels, with a Doppler frequency of 50 Hz and Sampling frequency of 100000 Hz


<img src="/Release/Error%20Correction%20Codes/Jakes_10Hz.png" width="49%">

The plot above shows the BER vs Eb/N0 plot for an correlated Rayleigh channel with a rate of 0.5, Doppler frequency of 50Hz, codeword length of 60 and Sampling frequency of 100000 Hz.



# Updated BER vs Eb/N0 peformance plots using LDPC Codes
The "Final" folder contains codes for using any parity check matrix from **https://www.inference.org.uk/mackay/codes/**, with any rate and codeword length, to examine the performance of LDPC coding on BPSK transmission using a variety of channels (AWGN, Rayleigh, Jakes). 

The LDPC Parity Check Matrix used is a systematic 408*204 matrix with rate 1/2 (204.33.484 (N=204,K=102,M=102,R=0.5)). 

The averaging has been performed for upto 10000000 bits. Reference plots are also given below. 

<img src="/Release/Error%20Correction%20Codes/AWGN.png" width="49%"> <img src="/Release/Error%20Correction%20Codes/AWGN_Mackay.png" width="41%">

The plot above on the left shows the BER vs Eb/N0 plot for an AWGN channel with a rate of 0.5 and codeword length of 408. The one on the right shows the expected performance as given on **https://www.inference.org.uk/mackay/codes/**. 

<img src="/Release/Error%20Correction%20Codes/Rayleigh.png" width="49%"> 

The plot above shows the BER vs Eb/N0 plot for a Rayleigh channel with a rate of 0.5 and codeword length of 408.

<img src="/Release/Error%20Correction%20Codes/Jakes100.png" width="41%">

The plot above on the left shows the BER vs Eb/N0 plot for a Jakes channel with a rate of 0.5 and codeword length of 408, with Doppler of 100 Hz and Fs = 10^6 Hz. 













 



  



