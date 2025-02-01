function [hw] = Jakes_filter(f_max, Ts, N)
% FIR channel shaping filter with Jakes Doppler spectrum
% S(f) = 1 / [sqrt(1 - (f/f_max)^2) * pi * f_max]
% Use impulse response obtained from spectral factorization
% Input parameters are
% f_max = maximum Doppler frequency in Hz
% Ts = sampling interval in seconds
% N = desired filter order
% Returns the windowed impulse response of the FIR filter

L = N / 2;
n = 1:L;
J_pos = besselj(0.25, 2 * pi * f_max * n * Ts) ./ (n.^0.25);
J_neg = fliplr(J_pos);
J_0 = 1.468813 * (f_max * Ts)^0.25;
J = [J_neg J_0 J_pos]; % Jakes filter
% Shaping filter smoothed using Hamming window
n = 0:N;
hamm = 0.54 - 0.46 * cos(2 * pi * n / N); % Hamming window
hw = J .* hamm; % Multiply Jakes filter with the Hamming window
hw = hw ./ sqrt(sum(abs(hw).^2)); % Normalized impulse response
end


%-- Colored noise generating with Jakes PSD---
fd=10; %maximum Doppler frequency
Fs=100;%sampling frequency
N=512;%order of Jakes filter
Ts=1/Fs;%sampling interval
[h] = Jakes_filter(fd,Ts,512);%design Jakes filter

%--Exciting the filter with white noise---
x=randn(1,10000); %Generate white noise of 10,000 samples
y = conv(x,h,'valid');%convolve with Jakes filter

%--Plot time domain and freq domain view of colored noise y
%2 Random Variables - Simulating Probabilistic Systems
subplot(1,2,1); plot(log10(abs(y))); axis tight; %Noise envelope
xlabel('n'); ylabel('y[n]');title('Noise samples')



