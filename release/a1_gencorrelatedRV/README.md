# Jakes Channel Generator

The objective of the first part of the project is to be able to generate channel instances which will be used to evaluate and train the models. The Jakes channel model is chosen due to its greater generalisation as opposed to rayleigh or rician fading models. It takes into account the multipath fading effects and also allows one to simulate different levels of doppler frequencies on the channel</br>

The jakes channel is generated in a sum of sinusoids method as outlined in <https://www.ee.iitm.ac.in/~giri/pdfs/EE5141/Jakes-Simulation.pdf>. This has been adapted from the Rapport textbook and can be easily implemented to generate a jakes channel.

## Preliminaries
The generation of jakes channel was preceded by a number of different steps which helped me understand random variables and gain an understading of how to generate them. The things done are:

1. Generate correlated random variables(count = 2) using the method outlined in the paper <https://www.sciencedirect.com/science/article/pii/S0266892017302539?ref=pdf_download&fr=RR-2&rr=887b3613691f936b>. The correlation plots were observed and helped gain intuition. These programs are found in **/programs/june9**.
1. A jakes channel was attempted to be generated in real time. This was done by generating temporally correlated random variables using ARMA processes and fed into a filter whose shape was identical to the doppler spectrum as seen in the Jakes model. The spatial correlation between channels(for instance, the antennas on a mobile phone will all have slightly different but correlated channel gains) were modelled using the matrix approach that was part of the paper in the previous part. In essence, two independent RVs can be correlated using a matrix multiplication. All of this is found in **/programs/june20**.
1. Two other models(Clarke and Rayleigh) are also simulated to see how they differ from Jakes and to gain an intuition for channel models. A number of problems were faced in the Jakes generation at this time(**detail required**). All of this is found in **/programs/july2**.
1. The method of generating the jakes spectrum(ARMA->Filter->Matrix correlation) seemed to not work as expected as the correlation plots do not look right. All of this is in **/programs/july7**. A switch is made from the filter based approach to the sum of sinusoids method.

## Final Program
The final program using the SOS method as mentioned above. The correlation is created using the same matrix multiplication methods as mentioned. This program has support for:
1. Both real and complex channels can be generated
1. Doppler frequency and sampling rate can be chosen as required. (**should I provide extra details?**).
1. The program can generalise and generate N correlated random variables whose correlation matrix is determined(**how is this matrix chosen?**). For testing, only one channel is modelled and correlation is not relevant.
1. Additional programs such as one that plots the correlation and the channel itself is also part of the **/final** directory. This will form the final generator for the project. **Importantly, the correlation and behavior of the channel match what is expected.**

## Plots
### Magnitude and autocorrelation for a real channel
![plot](/release/a1_gencorrelatedRV/final/figure1_real_channel_magnitude_bessel.png)

The above plot shows the generated jakes channel using the sum of sinusoids method. The plot on the right is used to show the autocorrelation function's similarity to the corresponding zeroth Bessel function.

### Magnitude and autocorrelation for complex channel

A plot very similart to the top is obtained if the real and imaginary parts of the complex channel's correlations are plotted. **How should the complex number itself be handled with the real and imaginary parts together?**

## Pending Work
The current code can generate a fixed number of samples for a given data point and dump it into a text file. The following improvements may be needed:
1. Make the generation streaming. The channel should be generated until the program is terminated.
1. The channel values(which are complex numbers) might need to be stored better instead of being stored in a txt file(maybe use it directly in the estimation and model training)
1. The code can be optimized to run faster and require less RAM. Ways should be looked up to figure how this could be done.