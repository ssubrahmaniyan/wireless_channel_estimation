This folder contains the full suite for the initial generation of samples for the prediction project.

The file jakes_totxt.py can be tweaked to choose the number of different channel random variables needed and the covariance and mean matrices can be customised as per needs(randomised for the generation part). The number of data points needed for the audio signal can be found out from the shape of the audio to be processed and fed into the generator file. It outputs a channel.txt file which contains the channel values thus generated using the jakes model.

The file fullprocess.py reads the audio signal and scales them with the channel values to generate the new signals. AWGN is also added and the cleaning of such a noise is shown using the LS separation method. 

In essence, this module with appropritate tweaks can be run many times on an audio signal to generate the input and output pairs that will be required for the sample generation for the prediction parts of the project.

Extension to do: All of this can be added into one script with appropriate arguments to process fully. For instance, the audio file, the number of different channel samples and noise ratio can be input in the script which then handles all the other aspects on its own. This is left to be done.
