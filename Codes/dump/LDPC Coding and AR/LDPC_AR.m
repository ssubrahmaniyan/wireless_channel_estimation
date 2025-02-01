% Main script for comparing all approaches
clear all;
close all;

% Common Parameters
EbNo_dB = -1:0.5:20;
N = 10000000;  % Total channel samples
max_frames = 100;
code_rate = 1/2;

% Load Jakes channel
channel_values = complex(zeros(N, 1));
fid = fopen('channels.txt', 'r');
for i = 1:N
    line = fgetl(fid);
    if ischar(line)
        values = jsondecode(line);
        channel_values(i) = complex(values(1), values(2));
    end
end
fclose(fid);

% Parameters
block_length = 324;  % Changed from 648 to 324 (uncoded block length)
rate = 1/2; % Code rate
max_iter = 50; % Maximum number of decoding iterations
var_order = 25;
pilot_block_size = block_length;  % Using block_length directly
num_initial_pilots = 2*var_order + 2;
mse_threshold = 0.1;

% Initialize LDPC code
ldpc = LDPCCode(2*block_length, block_length);  % Changed to use 2*block_length for N
ldpc.load_wifi_ldpc(2*block_length, rate);

% Initialize BER arrays
BER_coded = zeros(size(EbNo_dB));
BER_uncoded1 = zeros(size(EbNo_dB));
BER_uncoded2 = zeros(size(EbNo_dB));
BER_ar_uncoded = zeros(size(EbNo_dB));
BER_ar_coded = zeros(size(EbNo_dB));
retrans_freq_uncoded = zeros(size(EbNo_dB));
retrans_freq_coded = zeros(size(EbNo_dB));

% Main simulation loop
for i = 1:length(EbNo_dB)
    fprintf('Simulating Eb/No = %.1f dB\n', EbNo_dB(i));
    
    % Standard Jakes simulations (from LDPC_Jakes.m)
    [BER_coded(i), BER_uncoded1(i), BER_uncoded2(i)] = simulateJakes(channel_values, ...
        EbNo_dB(i), max_frames, ldpc);
    
    % AR Model simulations
    [retrans_freq_uncoded(i), BER_ar_uncoded(i)] = simulateAR(channel_values, ...
        false, pilot_block_size, var_order, mse_threshold, num_initial_pilots, EbNo_dB(i), ldpc);
    
    [retrans_freq_coded(i), BER_ar_coded(i)] = simulateAR(channel_values, ...
        true, pilot_block_size, var_order, mse_threshold, num_initial_pilots, EbNo_dB(i), ldpc);
end

% Plotting
figure(1);
semilogy(EbNo_dB, BER_coded, 'k-o', 'DisplayName', 'Coded Jakes');
hold on;
semilogy(EbNo_dB, BER_uncoded1, 'b-s', 'DisplayName', 'Uncoded Jakes');
semilogy(EbNo_dB, BER_uncoded2, 'r-^', 'DisplayName', 'Uncoded Rayleigh');
semilogy(EbNo_dB, BER_ar_uncoded, 'm-v', 'DisplayName', 'AR Uncoded');
semilogy(EbNo_dB, BER_ar_coded, 'g-d', 'DisplayName', 'AR Coded');
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER Performance Comparison');
legend('Location', 'southwest');

figure(2);
plot(EbNo_dB, retrans_freq_uncoded, 'r-o', 'DisplayName', 'AR Uncoded');
hold on;
plot(EbNo_dB, retrans_freq_coded, 'b-s', 'DisplayName', 'AR Coded');
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Retransmission Frequency');
title('Pilot Retransmission Frequency');
legend('Location', 'northeast');

% Helper function for AR simulation
function [retrans_freq, ber] = simulateAR(channel_values, use_coding, block_size, ...
    var_order, mse_threshold, num_initial_pilots, EbNo_dB, ldpc)
    
    % Initialize variables
    current_idx = 1;
    total_bits = 0;
    wrong_bits = 0;
    pilot_requests = 0;
    
    % Get initial channel estimates
    initial_bits = num_initial_pilots * block_size;
    channel_estimates = getInitialEstimates(channel_values, current_idx, ...
        num_initial_pilots, block_size, EbNo_dB);
    current_idx = current_idx + initial_bits;
    pilot_requests = pilot_requests + 1;
    
    while current_idx + block_size <= length(channel_values)
        % Train VAR model
        channel_data = [real(channel_estimates), imag(channel_estimates)];
        mdl = varm(2, var_order);
        mdl = estimate(mdl, channel_data);
        
        % Predict next channel value
        [Y_pred, ~] = forecast(mdl, 1, channel_data);
        htilde = complex(Y_pred(1), Y_pred(2));
        
        if use_coding
            [wrong_bits_block, hnought] = processCodedBlock(channel_values, current_idx, ...
                htilde, block_size, EbNo_dB, ldpc);
        else
            [wrong_bits_block, hnought] = processUncodedBlock(channel_values, current_idx, ...
                htilde, block_size, EbNo_dB);
        end
        
        % Check MSE
        mse = abs(htilde - hnought)^2;
        
        if mse > mse_threshold
            % Retransmit pilots
            channel_estimates = getInitialEstimates(channel_values, current_idx, ...
                num_initial_pilots, block_size, EbNo_dB);
            current_idx = current_idx + initial_bits;
            pilot_requests = pilot_requests + 1;
        else
            % Update and continue
            channel_estimates = [channel_estimates(2:end); hnought];
            wrong_bits = wrong_bits + wrong_bits_block;
            total_bits = total_bits + block_size;
            current_idx = current_idx + block_size;
        end
    end
    
    retrans_freq = pilot_requests * initial_bits / total_bits;
    ber = wrong_bits / total_bits;
end

function channel_estimates = getInitialEstimates(channel_values, current_idx, num_pilots, block_size, EbNo_dB)
    % Initialize channel estimates array
    channel_estimates = zeros(num_pilots, 1);
    
    % Convert Eb/No to SNR
    SNR = 10^(EbNo_dB/10);
    
    % For each pilot block
    for i = 1:num_pilots
        % Generate random pilot bits
        pilot_bits = randi([0 1], block_size, 1);
        symbols = 1 - 2*pilot_bits;
        
        % Get actual channel values for this block
        block_start = current_idx + (i-1)*block_size;
        block_end = block_start + block_size - 1;
        actual_channel = channel_values(block_start:block_end);
        
        % Add noise
        noise = sqrt(1/(2*SNR)) * (randn(block_size,1) + 1j*randn(block_size,1));
        received = actual_channel.*symbols + noise;
        
        % Estimate channel using least squares
        channel_estimates(i) = mean(received./symbols);
    end
end

function [wrong_bits, hnought] = processCodedBlock(channel_values, current_idx, htilde, block_size, EbNo_dB, ldpc)
    % Parameters
    code_rate = 1/2;
    snr_db = EbNo_dB + 10*log10(code_rate);
    SNR = 10^(snr_db/10);
    noise_var = 1/(2*SNR);
    
    % Generate random information bits
    info_bits = randi([0 1], block_size, 1);  % Now block_size is 324
    
    % Encode using LDPC
    coded_bits = ldpc.encode_bits(info_bits);  % Will produce 648 bits
    symbols = 1 - 2*coded_bits;
    
    % Get actual channel values for full coded block
    actual_channel = channel_values(current_idx:current_idx+2*block_size-1);  % Changed to 2*block_size
    
    % Add noise
    noise = sqrt(noise_var/2) * (randn(2*block_size,1) + 1j*randn(2*block_size,1));  % Changed size
    received = actual_channel.*symbols + noise;
    
    % Compute LLRs using predicted channel
    llr = 2*real(conj(htilde).*received)./(abs(htilde)^2 * noise_var);
    
    % Decode using LDPC
    decoded_bits = ldpc.decode_llr(llr, 50, true);
    
    % Count errors (only in information bits)
    wrong_bits = sum(decoded_bits(1:block_size) ~= info_bits);
    
    % Re-estimate channel using decoded bits
    hnought = mean(received./(1-2*coded_bits));
end

function [wrong_bits, hnought] = processUncodedBlock(channel_values, current_idx, htilde, block_size, EbNo_dB)
    % Convert Eb/No to SNR
    SNR = 10^(EbNo_dB/10);
    noise_var = 1/(2*SNR);
    
    % Generate random information bits
    info_bits = randi([0 1], block_size, 1);
    symbols = 1 - 2*info_bits;
    
    % Get actual channel values
    actual_channel = channel_values(current_idx:current_idx+block_size-1);
    
    % Add noise
    noise = sqrt(noise_var/2) * (randn(block_size,1) + 1j*randn(block_size,1));
    received = actual_channel.*symbols + noise;
    
    % Demodulate using predicted channel
    decisions = real(received.*conj(htilde))./(abs(htilde)^2) < 0;
    
    % Count errors
    wrong_bits = sum(decisions ~= info_bits);
    
    % Re-estimate channel using decisions
    hnought = mean(received./(1-2*decisions));
end

function [ber_coded, ber_uncoded1, ber_uncoded2] = simulateJakes(channel_values, EbNo_dB, max_frames, ldpc)
    % Initialize counters
    bit_errors_coded = 0;
    total_bits_coded = 0;
    bit_errors_uncoded1 = 0;
    total_bits_uncoded1 = 0;
    bit_errors_uncoded2 = 0;
    total_bits_uncoded2 = 0;
    
    % Parameters
    block_length = ldpc.K;  % This is now 324
    rate = 1/2;
    current_index = 1;
    
    % Convert Eb/No to SNR
    snr_db_coded = EbNo_dB + 10*log10(rate);
    SNR_coded = 10^(snr_db_coded/10);
    noise_var_coded = 1/(2*SNR_coded);
    
    SNR_uncoded = 10^(EbNo_dB/10);
    noise_var_uncoded = 1/(2*SNR_uncoded);
    
    for frame = 1:max_frames
        if current_index + 2*block_length > length(channel_values)
            break;
        end
        
        % Coded transmission
        info_bits = randi([0 1], ldpc.K, 1);
        codeword = ldpc.encode_bits(info_bits);
        symbols = 1 - 2*codeword;
        
        channel = channel_values(current_index:current_index + ldpc.N - 1);
        noise = sqrt(noise_var_coded/2) * (randn(ldpc.N,1) + 1j*randn(ldpc.N,1));
        received = channel.*symbols + noise;
        
        llr = 2*real(conj(channel).*received)./(abs(channel).^2 * noise_var_coded);
        decoded = ldpc.decode_llr(llr, 50, true);
        
        bit_errors_coded = bit_errors_coded + sum(decoded(1:ldpc.K) ~= info_bits);
        total_bits_coded = total_bits_coded + ldpc.K;
        
        current_index = current_index + ldpc.N;
        
        % Uncoded transmission 1 (with Jakes channel)
        info_bits = randi([0 1], block_length, 1);
        symbols = 1 - 2*info_bits;
        
        channel = channel_values(current_index:current_index + block_length - 1);
        noise = sqrt(noise_var_uncoded/2) * (randn(block_length,1) + 1j*randn(block_length,1));
        received = channel.*symbols + noise;
        
        decisions = real(received./channel) < 0;
        
        bit_errors_uncoded1 = bit_errors_uncoded1 + sum(decisions ~= info_bits);
        total_bits_uncoded1 = total_bits_uncoded1 + block_length;
        
        current_index = current_index + block_length;
        
        % Uncoded transmission 2 (with Rayleigh channel)
        info_bits = randi([0 1], block_length, 1);
        symbols = 1 - 2*info_bits;
        
        channel = sqrt(0.5) * (randn(block_length,1) + 1j*randn(block_length,1));
        noise = sqrt(noise_var_uncoded/2) * (randn(block_length,1) + 1j*randn(block_length,1));
        received = channel.*symbols + noise;
        
        decisions = real(received./channel) < 0;
        
        bit_errors_uncoded2 = bit_errors_uncoded2 + sum(decisions ~= info_bits);
        total_bits_uncoded2 = total_bits_uncoded2 + block_length;
    end
    
    ber_coded = bit_errors_coded / total_bits_coded;
    ber_uncoded1 = bit_errors_uncoded1 / total_bits_uncoded1;
    ber_uncoded2 = bit_errors_uncoded2 / total_bits_uncoded2;
end

