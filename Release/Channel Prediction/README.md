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

I also need to edit the code to train the VAR model on all the estimated channel values, rather than the channel values that are correctly (i.e. MSE is within the threshold) predicted, and observe the results. 

### Plots
![hi](/Release/Channel%20Prediction/figure1_real_channel_magnitude_bessel.png)
The plot above shows the behaviour of the 

![hi](/Release/Channel%20Prediction/figure1_real_channel_magnitude_bessel.png)
