# Autoregressive Prediction of Channels

This part of the project aims to predict the channel values using autoregressive models. 

## Simple Autoregressive Model

The first attempt is to use a simple autoregressive model to predict the channel. The idea is executed in steps as follows:
1. Send a few packets of pilot bits and use an estimation method to gain estimates for the channels.
1. Train a VAR model(chosen because the channel is complex) to learn the channels and make a prediction for the next channel value.
1. Use the sign of this channel value(because the scheme used is BPSK) and find the data that was sent. Use this data and an estimation method to estimate the true channel value(note that this is different from the predicted. Estimated and predicted are different)/
1. Check if this prediction is good enough by computing the MSE. If the MSE exceeds a threshold, resend a few pilot packets and restart from step 1. Otherwise, append the new prediction to the model arguments, refit and make another prediction and continue the process.


### Results
This method was performed using a bash script to run on different MSE thresholds for different values and the plot of the pilot ratio against the normalised doppler frequency(Fd/Fs) is shown below. **As expected, the values saturate as the doppler increases and a lower MSE threshold requires a larger number of pilots and hence a larger pilot ratio.**
The data is stored in the file parameter_sweep_data.txt

![plot](/release/b1_prediction/b11_autoregression/simpleVAR/ratio_vs_doppler.png)
For this plot, the script was run with Fs = 100000 and number of points = 500000 per simulation.
The x axis is the ratio Fd/Fs as asked for.

![plot](/release/b1_prediction/b11_autoregression/simpleVAR/ratio_vs_SNR_plot.png)

Same plot as the above, but with a variation in the SNR. The number of points is only 100000 here though. Note that the ratio is not plotted on the x axis and the true doppler values are printed.


### Questions/Extension work:
1. Can we use a model other than least squares to estimate values? (The suggested was a robust linear estimator).
1. Tune the parameters such as pilot size, and var_order to see performance differences. The coherence time approximation tells us that the pilot packet size must be less than k * Fs/Fd. Should this be used appropriately depending on the doppler that we aim to support in our system?
**When this is made use of, the plot seems to be a little different than others, this is located in the tempimages directory**

1. How can I make the VAR model work better than it already does? (Apart from parameter tuning, is there anything else to do?)
1. Everytime the MSE threshold is crossed, how many pilots should be sent to recalibrate the model? Should we do it as: send one, refit, check. If better, continue otherwise send another refit and check?

## Organisational Task

1. Should I instead divide the files into images, data and scripts in different folders and store them?

## Team meet questions:

1. Here is another approach I thought about: How about I try running the simulations for different packet/data ratios and compare the percentages of the packets which can be predicted with a specified mse accuracy(for the header bits) and a specified ber accuracy(for the data bits)
