# Autoregressive Prediction of Channels

This part of the project aims to predict the channel values using autoregressive models. 

**The simple autoregressive model, which itself has is data driven(I have named the other model as data driven because I simply felt like doing so) is the one with coherent results and hence would be considered final model for further progress**

## Simple Autoregressive Model

**All of the programs, plots and data relevant to this model is in the folder /simpleVAR in the same directory**

The first attempt is to use a simple autoregressive model(with a data driven approach) to predict the channel. The idea is executed in steps as follows:
1. Send a few packets of pilot bits and use an estimation method to gain estimates for the channels.
1. Train a VAR model(chosen because the channel is complex) to learn the channels and make a prediction for the next channel value.
1. Use the sign of this channel value(because the scheme used is BPSK) and find the data that was sent. Use this data and an estimation method to estimate the true channel value(note that this is different from the predicted. Estimated and predicted are different)/
1. Check if this prediction is good enough by computing the MSE. If the MSE exceeds a threshold, resend a few pilot packets(varying the number of such packets will vary the results) and restart from step 1. Otherwise, append the new prediction to the model arguments, refit and make another prediction and continue the process.


### Results
This method was performed using a bash script to run on different MSE thresholds for different values and the plot of the pilot ratio against the normalised doppler frequency(Fd/Fs) is shown below. **As expected, the values saturate as the doppler increases and a lower MSE threshold requires a larger number of pilots and hence a larger pilot ratio.**
The data is stored in the file parameter_sweep_data.txt

![plot](/release/b1_prediction/b11_autoregression/simpleVAR/ratio_vs_doppler.png)
For this plot, the script was run with Fs = 100000 and number of points = 500000 per simulation.
The x axis is the ratio Fd/Fs as asked for.

![plot](/release/b1_prediction/b11_autoregression/simpleVAR/ratio_vs_SNR_plot.png)

Same plot as the above, but with a variation in the SNR. The number of points is only 100000 here though. Note that the ratio is not plotted on the x axis and the true doppler values are printed.

Both of these plots are made such that whenever the MSE threshold is exceeded, the numbers of pilots sent to recalibrate the VAR model is var_order + 2, which essentially resets the model altogether.

### Questions/Extension work:
1. Can we use a model other than least squares to estimate values? (The suggested was a robust linear estimator).
1. Tune the parameters such as pilot size, and var_order to see performance differences. The coherence time approximation tells us that the pilot packet size must be less than k * Fs/Fd. Should this be used appropriately depending on the doppler that we aim to support in our system?
~~**When this is made use of, the plot seems to be a little different than others, this is located in the tempimages directory**~~

1. How can I make the VAR model work better than it already does? (Apart from parameter tuning, is there anything else to do?)
1. Everytime the MSE threshold is crossed, how many pilots should be sent to recalibrate the model? Should we do it as: send one, refit, check. If better, continue otherwise send another refit and check?

## Organisational Task

~~1. Should I instead divide the files into images, data and scripts in different folders and store them?~~

## Team meet questions:

~~1. Here is another approach I thought about: How about I try running the simulations for different packet/data ratios and compare the percentages of the packets which can be predicted with a specified mse accuracy(for the header bits) and a specified ber accuracy(for the data bits)~~


## Modified Autoregressive Model

**All of the programs, plots and data related to this model are in the directory /datadrivenVAR**

This model is very similar to the previous implementation except in that at an instant, either a pilot packet may be requested or a whole conventional packet(consisting of some pilot and some payload bits) may be requested.

The process implemented in this part is:
1. Send a pilot packets to tune the VAR model.
1. Make a prediction for the channel value of the next quanta. I have called the size of a pilot packet a quanta because VAR models require the data to be sampled in uniform time - the samples should be of the same size and hence all the inputs must be equall separated.
1. Send a conventional packet and estimate the channel value for the pilot part of the packet. Add this estimate to the training data for the VAR model.
1. Refit the model and make another prediction for the next data quanta. Use the prediction to demodulate the data quanta, gain a refined estimate for the channel value for that quanta and retrain and repeat until all the data quanta are done.
1. When the whole packet is processed, compute the MSE between the predicted and estimated h values for all of the quantas in the packet. If this is greater than the threshold, send a bunch of pilots to refit the model and repeat the same process.
1. The BER values are computed for the data packets as done usually.

By running this simulation for different values of MSE thresholds, pilot size and quanta of data packets involved, a spot at which the BER is minimum and the pilot ratio is minimum can be attained. This point can be found out by running numerous iterations of the same program - **yet to be done properly**

![plot](/release/b1_prediction/b11_autoregression/datadrivenVAR/datadrivenVARsweep.png)

The configurations on the x axis can be understood by looking at the results.txt file. The parameters used in the program are: <br>
Fd = 100
Fs = 100000
snr_db = 10
var_order = 25
initial_pilots = var_order + 2
</br>


This plot is generated by running the sweep_datadrivenVAR.py file in the same directory. The resulting data is stored in the results.txt file. Note that the behavior is as expected:

1. As the MSE threshold increases for a fixed pilot size and a fixed quanta, the BER remains aproximately the same and the pilot ratio decreases.
1. As the quanta increases for a fixed MSE threshold and pilot size, the BER increases and the pilot ratio decreases as expected.

## Things left to do
1. The plot for a variable packet sizes must be redone.
1. The modified data driven model's plot must be made clearer to make sense.
1. Fix how the refined estimates are made on the interpolation and non interpolation models.
1. Add the interpolation method alongside the VAR method for predicting and estimating the channels.

# Effect of positions of pilots in channel estimation

## Simulation of simple channel transmission
To demonstrate the working of the channel generation as expected, the generated Jakes channel values are used to perform a simple wireless communication method: data is modulated, transmit, demodulated and the errors are counted. The assumption here is that the channel values are known for each of the sent data packets. The corresponding plots are generated using the **theoretical.py** file.

![plot](/release/b1_prediction/b12_pilot_positioning/1Hz_theoretical.png)
This plot is created with 100000 complex channel points, Fd = 1 Hz and all other parameters default to the **jakeschannel.py** file.

![plot](/release/b1_prediction/b12_pilot_positioning/100Hz_theoretical.png)
This plot is created with 100000 complex channel points, Fd = 100 Hz and all other parameters are default in the file.

*The non smoothness of the plots are likely because of insufficient averaging in transmission. Using more points would give smoother plots as would be expected. Also, ignore the title in the images, the plots represent BER for BPSK in Jakes fading channels of different doppler frequency support alongside the theoretical plots for Rayleigh fading(corresponding to a doppler frequency of infinity and the AWGN channel corresponding to a doppler frequency of 0*

## Simulation of simple transmission with a data-driven approach

Here, the simulation is made packetwise - a pilot and payload are sent, the channel value is estimated from the pilot and used to demodulate the payload. Using the demodulated payload, we re-figure the channel values for the payload part alone, and again demodulate the received payload data. 

The files for this part of the implementation can be found in b12_pilot_positioning/nointerpolation

### Channel with a Doppler of 1 Hz

For the method outlined above, the following plot is obtained. The data corresponding to the same plot are found in 1Hz_datadriven_nointerpolation.txt.

![plot](/release/b1_prediction/b12_pilot_positioning/nointerpolation/1Hz_datadriven_nointerpolation.png)

The parameters with which the simulation was performed are a total of 100000 data points, packet_size = pilot_size = 5 bits. 

*The plot shows the pilot ratio used to attain the corresponding BER for a given MSE threshold, alongside the average MSE between the prediction and assumption of channel values(essentially, when using a data driven model, the error is between the channel values that are found out using the demodulated data and the one that is assumed to be constant from the previous packet).*

### Channel with Doppler of 100 Hz

For all the conditions identical to the one for 1Hz, except with a Doppler of 100Hz, the plot obtained is:

![plot](/release/b1_prediction/b12_pilot_positioning/nointerpolation/100Hz_datadriven_nointerpolation.png))

The corresponding data is found in the csv file with the same name.

### Observations
1. It is noted that as the packet and pilot sizes are increased, the BER curves flatten out at higher and higher BER. This is probably because the assumption of a constant channel value from the previous packet becomes a bad assumption.
1. The MSE, pilot ratio and BER curves all fall with an increasing SNR(same as Eb/N0 for this case), which is the expected behavior.
