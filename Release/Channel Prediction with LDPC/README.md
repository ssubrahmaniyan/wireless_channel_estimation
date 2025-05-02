# Channel Prediction with LDPC

Below is the plot for benchmarking performance of LDPC encoding and decoding on an AWGN channel
<img src="/Release/Channel%20Prediction%20with%20LDPC/LDPC_AR.png" width="48%"> <img src="/Release/Channel%20Prediction%20with%20LDPC/LDPC_AR_Sionna.png" width="51%"> 
The plot on the left is the one obtained by my simulation, and the one on the right is the refernece plot which can be obtained at **https://nvlabs.github.io/sionna/examples/Sionna_tutorial_part1.html**

Below is the plot for performance of an Autoregression-based channel prediction approach which is data-driven (MSE between estimated and predicted channels is used to determine if a retransmission of pilots is needed)
 <img src="/Release/Channel%20Prediction%20with%20LDPC/Uncoded_RetransmissionFreq_BER.png" width="80%">

 The code for reproducing the plots obtained above, as well as accessing the LDPC encoding and decoding functions for Autoregression approaches that use LDPC for FEC, can be found in the notebook LDPC_AR.ipynb

## Some observations after making corrections to earlier VAR attempts
Corrections made : 
* I am now including the BER of data packets that satisfy the pilot retransmission condition as well (I was previously ignoring such data packets). This leads to a minor increase in the BER for the iteration.
* An error in the Autoregression code so far is that we considered the data bits to count for errors (which is correct), but we divided by the total (data + pilot) bits transmitted to compute the BER. After fixing the error, the BER values are much higher, though the benefits of LDPC are still evident.


For the below parameters:
* Fs = 100000  # Sampling frequency
* Fd = 10    # Doppler frequency
* mse_threshold = 0.1  # Threshold for MSE
* N = 100000   # Total samples
* ebno_db = 10     # Signal-to-noise ratio in dB (interpreted as Eb/N0 for coded system)
* packet_size = 50  # Bits per packet
* pilot_size = 50  # Pilot bits
* var_order = 25  # VAR model order
* initial_pilots = 27  # Initial pilots

For investigating the effect of the block length on the LDPC performance, here is the plot for the variation of the channel for Fd = 100 Hz and Fd = 10 Hz with Fs = 100000 Hz

<img src="/Release/Channel%20Prediction%20with%20LDPC/100Hz_Channel.png" width="40%"> <img src="/Release/Channel%20Prediction%20with%20LDPC/10Hz_Channel.png" width="40%"> 

The outputs from 3 different approaches are:
| Algorithm | Retransmission Frequency | BER | Total Data Bits Transmitted (out of 100000) |
|-----------|--------------------------|-----|-------------------------------------------|
| Data Driven (decisions based on MSE) Channel Prediction with Autoregression without LDPC | 0.014496 | 0.01555 | 98600 |
| Data driven (decisions based on LLR) Channel Prediction with Autoregression with  LDPC | 0.253173828125 | 0.057342529296875 | 42750 |
| Data driven (decisions based on LLR with CRC checks) Channel Prediction with Autoregression with LDPC | 0.19091796875 | 0.0243682861328125 | 34750 |


Below is also a plot of the correlation between LLRs and Bit Errors in a LDPC codeword length of 60 and code rate of 0.5 (this is the first transmitted codeword)

<img src="/Release/Channel%20Prediction%20with%20LDPC/LLR_Error_Plot.png" width="80%">

************************************************************************************************************************************************************************************************
For the below parameters:
* Fs = 1000000  # Sampling frequency
* Fd = 10    # Doppler frequency
* mse_threshold = 0.1  # Threshold for MSE
* N = 1000000   # Total samples
* ebno_db = 10     # Signal-to-noise ratio in dB (interpreted as Eb/N0 for coded system)
* packet_size = 500  # Bits per packet
* pilot_size = 500  # Pilot bits #around packet size, ar 15, Fs 1000000, Fd 20, 10^6 * 100, 
* var_order = 15  # VAR model order
* initial_pilots = 27  # Initial pilots

The following changes have been implemented for the 3 simulations reported in the table below:
1. The Sum of Sinusoids method is now being used for generating the channel
2. The simulations are averaged 10 times and run for N = 100000 samples
3. The LLR cutoff is set to 12 and the cutoff frequency is set to 0.4
4. The VAR order has been reduced and the pilot size is now comparable to the size of the packet for which the channel is predicted. 


The outputs from 3 different approaches are:
| Algorithm | Retransmission Frequency | BER | Total Data Bits Transmitted (out of 100000) |
|-----------|--------------------------|-----|-------------------------------------------|
| Data Driven (decisions based on MSE) Channel Prediction with Autoregression without LDPC | 0.4695 | 0.0002295 | 61000 |
| Data driven (decisions based on LLR) Channel Prediction with Autoregression with LDPC | 0.23889 | 0.0000365 | 43000 |
| Data driven (decisions based on LLR with CRC checks) Channel Prediction with Autoregression with LDPC | 0.2668 | 0.000021875 | 40950 | 


Attached below are plots of the LLR magnitudes based on which a decision was made, and the BER vs Doppler frequency plot for the LDPC (without CRC) coded AR transmission. 

<img src="/Release/Channel%20Prediction%20with%20LDPC/LLR_Plot.png" width="45%"> <img src="/Release/Channel%20Prediction%20with%20LDPC/BER_vs_Doppler_P1.png" width="45%"> 

************************************************************************************************************************************************************************************************

 Block Error Rates (BLER) are now used in the transmission simulation. Each LDPC block with CRC bits sent is considered as a block, and CRC is used to verify if the block has been received correctly or not. The number of wrong blocks received gives us the BLER. 

For the below parameters:
* Fs = 1000000  Hz # Sampling frequency
* Fd = 10 to 100 Hz # Doppler frequency
* N = 1000000   # Total samples
* ebno_db = 10     # Signal-to-noise ratio in dB (interpreted as Eb/N0 for coded system)
* packet_size = 500  # Bits per packet
* pilot_size = 500  # Pilot bits #around packet size, ar 15, Fs 1000000, Fd 20, 10^6 * 100, 
* var_order = 15  # VAR model order
* initial_pilots = 27  # Initial pilots

The outputs from 2 different approaches for Fd = 100 Hz are:
| Algorithm | Retransmission Frequency | BER | Total Data Bits Transmitted (out of 100000) |
|-----------|--------------------------|-----|-------------------------------------------|
| Transmission with perfect CSI and LDPC FEC | NA | 0.00086 | NA |
| Data driven (decisions based on LLR) Channel Prediction with Autoregression using LDPC | 0.2668 | 0.00124 | 45150 |

The outputs from 2 different approaches for Fd = 100 Hz are:
| Algorithm | Retransmission Frequency | BLER | Total Data Bits Transmitted (out of 100000) |
|-----------|--------------------------|-----|-------------------------------------------|
| Transmission with perfect CSI and LDPC FEC (BLER) | NA | 0.08 | NA |
| Data driven (decisions based on CRC checks) Channel Prediction with Autoregression using LDPC (BLER) | 0.2712 | 0.115 | 40500 | 

Attached below are plots of the BER vs Doppler Frequency plots for the BER and BLER simulation detailed above, respectively

<img src="/Release/Channel%20Prediction%20with%20LDPC/BER_plot_LDPC.png" width="45%"> <img src="/Release/Channel%20Prediction%20with%20LDPC/BLER_plot_LDPC.png" width="45%"> 
 
