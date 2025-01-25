# Channel Prediction using Autoregressive models
### Algorithm 
The process of sending header bits and estimating channels for the entire packet length leads to a reduction in the effective number of data bits that can be reliably transmitted. Hence, a method used to reduce the number of headers in the trasmission is to predict channel values into the future. 
This is an attempt to train a Vector Autoregressive (VAR) model with channel values and then predict future channel values, following which I evaluate the accuracy of the predicted channels and the reliability of the predicted data bits. 

The algorithm used is as follows:
* Begin the transmission by sending a sequence of header packets, and estimate a sequence of channel values using the Least Squares method
* Train a VAR model on the sequence of the channel values determined
* Predict a few channel values into the future, and use these as channel estimates for data packets being sent.
* Demodulate the data packets with predicted channels and obtain a prediction on the bits being transmitted.
* Get another estimate of the channel value by performing Least Squares estimation on the estimated bits from the previous step and the received bits.
* Compute the MSE between the predicted and estimated channel values.

The purpose behind computing the MSE values is to determine when to retransmit header packets and retrain the VAR model with a fresh and accurate set of channel values. We have examined an arbitrary range of MSE thresholds from 0.01 to 0.1, and when the MSE for any set of predicted channels exceeds this threshold, we request a retransmission of header packets, thereby ensuring accuracy in the predictions of the VAR model. 

An evaluation metric to evaluate the performance of the VAR model is to compute the retransmission frequency of header packets, which is simply the number of pilot packets transmitted divided by the total number of packets transmitted. This is necessary since the aim of using VAR is to reduce the number of header packets needed for reliable communication. 

I also need to edit the code to train the VAR model on all the estimated channel values, rather than the channel values that are correctly (i.e. MSE is within the threshold) predicted, and observe the results. 

### Plots
<img src="/Release/Channel%20Prediction/Retransmission_Freq_Vs_Doppler_VAR.png" width="49%">
The plot above shows the behaviour of the VAR model characterized by the Retransmission Frequency computed as above, plotted against the Doppler frequency. Note that as the Doppler frequency increases, the channel values become harder to predict and the need to retransmit haders increases. Note also that this is at a fixed threshold of MSE. 

<img src="/Release/Channel%20Prediction/Retransmission_Freq_Vs_Doppler_VAR_Thresholds.png" width="49%">
This plot examines the performance of the VAR model for a correlated Rayleigh channel with Jakes spectrum at different MSE thresholds from 0.01 to 0.1, and clearly the retransmissions increase at a given Doppler frequency for a lower MSE threshold. 

The notebook Channel_Estimation_Pipeline_outputs.ipynb contains code to reproduce these outputs from scratch, and I need to fix the last plot generated in the notebook relating to BER Vs SNR. 

# Channel Prediction with Autoregressive models and LDPC coding

LDPC coding will increase the accuracy of the decoding process. For reference, the gains observed for a correlated Rayleigh channel are given below:

The plot below descirbes the coding gains obtained for the following parameters: Fs = 100000 Hz, Fd = 100 Hz, N = 10000000, LDPC Block length = 648

<img src="/Release/Channel%20Prediction/Jakes_Rayleigh2.png" width="49%">

The plot below describes the coding gains obtained for the following parameters: Fs = 100000 Hz, Fd = 10 Hz, N = 10000000, LDPC Block length = 1296

<img src="/Release/Channel%20Prediction/Jakes_Rayleigh.png" width="49%">

The plot below descirbes the coding gains obtained for the following parameters: Fs = 100000 Hz, Fd = 1 Hz, N = 10000000, LDPC Block length = 648

<img src="/Release/Channel%20Prediction/Jakes_Rayleigh3.png" width="49%">

Approaches evaluated for predicting channel values:
- Running the simulation bit-by-bit as done in previous AR simulations
- Using smaller block lengths with custom parity check matrices - led to poor demodulation results for standard channels
- Using specific block lengths in pilot and data packet transmission, training the VAR on the least squares estimate of a sequence of channel values.

The plot below descirbes the coding gains obtained for the following parameters: Fs = 10000 Hz, Fd = 15.43 Hz, N = 10000000, LDPC Block length = 648

<img src="/Release/Channel%20Prediction/Jakes_AR.png" width="49%">

The plot below describes the coding gains obtained for the following parameters: Fs = 100000 Hz, Fd = 10 Hz, N = 10000000, LDPC Block length = 1296, as obtained before. This block length was chosen to ensure that several coherence times are spanned for a wide range of Doppler frequencies and also because of the quality of codes available at this length. 

<img src="/Release/Channel%20Prediction/Jakes_Rayleigh.png" width="49%">

The plot below describes the channel variation for the above mentioned parameters: Fs = 100000 Hz, Fd = 10 Hz, N = 10000000, LDPC Block length = 1296. It shows the accuracy of the VAR model in MATLAB in predicting the channel values given sufficient training data (order 25, 200 actual samples being trained on)

<img src="/Release/Channel%20Prediction/10000_10_10000000.png" width="49%">

Issues faced in using alist format codes:
1. Conversion from alist to MATLAB sparse parity check matrices is completed
2. Since my codes run for specific block lengths, I attempted to use ldpcEncode and ldpcDecode functions from MATLAB.
3. Performance on different combinations of Fs, Fd and parity check matrices was very poor, and on Mackay codes, an error was raised saying the dimensions of the codes are incorrect. 
4. Need to check the format of the H matrices needed for the MATLAB defined functions. 








