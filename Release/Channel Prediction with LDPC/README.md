# Channel Prediction with LDPC

Below is the plot for benchmarking performance of LDPC encoding and decoding on an AWGN channel
<img src="/Release/Channel%20Prediction%20with%20LDPC/LDPC_AR.png" width="49%"> <img src="/Release/Channel%20Prediction%20with%20LDPC/LDPC_AR_Sionna.png" width="49%"> 
The plot on the left is the one obtained by my simulation, and the one on the right is the refernece plot which can be obtained at **https://nvlabs.github.io/sionna/examples/Sionna_tutorial_part1.html**

Below is the plot for performance of an Autoregression-based channel prediction approach which is data-driven (MSE between estimated and predicted channels is used to determine if a retransmission of pilots is needed)
 <img src="![Uncoded_RetransmissionFreq_BER](https://github.com/user-attachments/assets/d31e91b4-a0ae-4a81-be9d-4e63b9678fdf)" width="60%">

 The code for reproducing the plots obtained above, as well as accessing the LDPC encoding and decoding functions for Autoregression approaches that use LDPC for FEC, can be found in the notebook LDPC_AR.ipynb
 

 
