# Wireless_Channel_Estimation

Workflow for SISO wireless communication

Generate a signal to be transmitted, with some bits in each packet as reference bits
Receive the modified packet at the receiver
Denoise the packet using least squares estimation, y = Ax, find an optimal x so MSE(y, Ax) is minimised
Use the reference bits and determine the channel for the packet (assumed to be time-invariant for the duration of the packet)
Use this channel to then infer the remaining of the original signal from the modified packet, iterate over whole signal
Proposed Alternate Workflow

Generate a signal to be transmitted, in the form of packets, but no header bits are labelled
Receive the modified packet at the receiver
Denoise the packet using least squares estimation as before
Given the sequence of prior channel values for prior packets received, use a sequence model (LSTM, GRU, Attention, GNN) to predict the current channel value
Use this channel to then infer the original packet from the modified packet, iterate over whole signal
The channel for metropolitan contexts adopts a Jakes' spectrum and Jakes' model of propagation, so the envelop of the channel will be a Rayleigh variable.
