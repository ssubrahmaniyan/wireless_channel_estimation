clear all;
close all;

% Paameters
block_length = 648; % Block length
rate = 1/2; % Code rate
max_iter = 50; % Maximum number of decoding iterations
EbNo_dB = -1:1:13; % Eb/No values in dB
max_frames = 100;

% Read channel values from file
N = 10000000;  % Total number of points
channel_values = complex(zeros(N, 1));  % Preallocate complex array
fid = fopen('channels.txt', 'r');
for i = 1:N
    line = fgetl(fid);
    if ischar(line)
        values = jsondecode(line);
        channel_values(i) = complex(values(1), values(2));
    end
end
fclose(fid);

% Initialize LDPC code
ldpc = LDPCCode(block_length, block_length * rate);
ldpc.load_wifi_ldpc(block_length, rate);

% Preallocate BER arrays
BER_coded = zeros(size(EbNo_dB));
BER_uncoded1 = zeros(size(EbNo_dB));
BER_uncoded2 = zeros(size(EbNo_dB));

% Convert Eb/No to SNR
snr_db_vec = EbNo_dB + 10*log10(rate);

% Track current position in channel values
current_index = 1;

% Loop over Eb/No values
for i = 1:length(EbNo_dB)
    fprintf('Simulating Eb/No = %.1f dB\n', EbNo_dB(i));
    SNR_linear = 10^(snr_db_vec(i)/10);
    noise_var = 1 / (2 * rate * SNR_linear);
    snr = 10^(snr_db_vec(i)/10);
    
    bit_errors_coded = 0;
    total_bits_coded = 0;
    bit_errors_uncoded1 = 0;
    total_bits_uncoded1 = 0;
    bit_errors_uncoded2 = 0;
    total_bits_uncoded2 = 0;
    frame = 0;
    
    while (frame < max_frames && current_index + 2*block_length <= length(channel_values))
        frame = frame + 1;
        
        %% Coded Transmission
        % Generate random information bits
        info_bits_coded = randi([0, 1], ldpc.K, 1);
        
        % Encode the bits
        codeword = ldpc.encode_bits(info_bits_coded);
        
        % BPSK modulation
        symbols_coded = 1 - 2 * codeword;
        
        % Get channel values for this frame
        channel_coded = channel_values(current_index:current_index + ldpc.N - 1);
        current_index = current_index + ldpc.N;
        
        % Generate complex noise
        noise_coded = sqrt(noise_var/2) * (randn(ldpc.N, 1) + 1j * randn(ldpc.N, 1));
        
        % Pass through channel with noise
        received_coded = channel_coded .* symbols_coded + noise_coded/sqrt(snr);
        
        % Compute LLRs
        llr_coded = real(2 * conj(channel_coded) .* received_coded ./ (abs(channel_coded).^2 * noise_var));
        
        % Decode using LDPC
        decoded_codeword = ldpc.decode_llr(llr_coded, max_iter, true);
        
        % Count bit errors for coded transmission
        bit_errors_coded = bit_errors_coded + sum(decoded_codeword(1:ldpc.K) ~= info_bits_coded);
        total_bits_coded = total_bits_coded + ldpc.K;
        
        %% Uncoded Transmission 1
        SNR_linear = 10^(EbNo_dB(i)/10);
        noise_var = 1 / (2 * SNR_linear);
        % Generate random information bits
        info_bits_uncoded = randi([0, 1], block_length, 1);
        
        % BPSK modulation
        symbols_uncoded = 1 - 2 * info_bits_uncoded;
        
        % Get channel values for this frame
        channel_uncoded = channel_values(current_index:current_index + block_length - 1);
        current_index = current_index + block_length;
        
        % Generate complex noise
        noise_uncoded = sqrt(noise_var/2) * (randn(block_length, 1) + 1j * randn(block_length, 1));
        
        % Pass through channel with noise
        received_uncoded = channel_uncoded .* symbols_uncoded + noise_uncoded/sqrt(SNR_linear);
        
        % Hard decision decoding
        decoded_uncoded = real(received_uncoded ./ channel_uncoded) < 0;
        
        % Count bit errors for uncoded transmission
        bit_errors_uncoded1 = bit_errors_uncoded1 + sum(decoded_uncoded ~= info_bits_uncoded);
        total_bits_uncoded1 = total_bits_uncoded1 + block_length;

        %%Uncoded Transmission 2
        % Generate random information bits
        info_bits_uncoded = randi([0, 1], block_length, 1);

        % BPSK modulation
        symbols_uncoded = 1 - 2 * info_bits_uncoded;

        % Generate Rayleigh channel and complex noise
        channel_uncoded = sqrt(0.5) * (randn(block_length, 1) + 1j * randn(block_length, 1));
        noise_uncoded = sqrt(1 / 2) * (randn(block_length, 1) + 1j * randn(block_length, 1));

        % Pass through channel with noise
        received_uncoded = channel_uncoded .* symbols_uncoded + (noise_uncoded/sqrt(SNR_linear));

        % Hard decision decoding using real part
        decoded_uncoded = real(received_uncoded ./ channel_uncoded) < 0;

        % Count bit errors for uncoded transmission
        bit_errors_uncoded2 = bit_errors_uncoded2 + sum(decoded_uncoded ~= info_bits_uncoded);
        total_bits_uncoded2 = total_bits_uncoded2 + block_length;
        
        % Print progress every 1000 frames
        if mod(frame, 1000) == 0
            fprintf('Frame %d: Coded BER = %.2e, Uncoded BER = %.2e\n', ...
                frame, bit_errors_coded/total_bits_coded, bit_errors_uncoded1/total_bits_uncoded1);
        end
    end
    
    % Reset channel index for next SNR point
    current_index = 1;
    
    % Calculate BER for the current Eb/No value
    BER_coded(i) = bit_errors_coded / total_bits_coded;
    BER_uncoded1(i) = bit_errors_uncoded1 / total_bits_uncoded1;
    BER_uncoded2(i) = bit_errors_uncoded2 / total_bits_uncoded2;
end

% Now initialize arrays for AR simulation
BER_ar_uncoded = zeros(size(EbNo_dB));
BER_ar_coded = zeros(size(EbNo_dB));
retrans_freq_uncoded = zeros(size(EbNo_dB));
retrans_freq_coded = zeros(size(EbNo_dB));

% AR parameters
var_order = 10;
num_initial_pilots = 4*var_order + 2;
mse_threshold = 0.1;

% Run AR simulations
for i = 1:length(EbNo_dB)
    fprintf('Simulating AR for Eb/No = %.1f dB\n', EbNo_dB(i));
    
    [retrans_freq_uncoded(i), BER_ar_uncoded(i)] = simulateAR(channel_values, ...
        false, block_length, var_order, mse_threshold, num_initial_pilots, EbNo_dB(i), ldpc);
    
    [retrans_freq_coded(i), BER_ar_coded(i)] = simulateAR(channel_values, ...
        true, block_length, var_order, mse_threshold, num_initial_pilots, EbNo_dB(i), ldpc);
end

% Plot all results
figure(1);
semilogy(EbNo_dB, BER_coded, '-o', 'LineWidth', 2, 'DisplayName', 'Coded');
hold on;
semilogy(EbNo_dB, BER_uncoded1, '-s', 'LineWidth', 2, 'DisplayName', 'Uncoded 1');
semilogy(EbNo_dB, BER_uncoded2, '-s', 'LineWidth', 2, 'DisplayName', 'Uncoded 2');
semilogy(EbNo_dB, BER_ar_uncoded, '-v', 'LineWidth', 2, 'DisplayName', 'AR Uncoded');
semilogy(EbNo_dB, BER_ar_coded, '-d', 'LineWidth', 2, 'DisplayName', 'AR Coded');
grid on;
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs Eb/No for Coded and Uncoded Transmission');
legend show;

figure(2);
plot(EbNo_dB, retrans_freq_uncoded, '-o', 'DisplayName', 'AR Uncoded');
hold on;
plot(EbNo_dB, retrans_freq_coded, '-s', 'DisplayName', 'AR Coded');
grid on;
xlabel('Eb/No (dB)');
ylabel('Retransmission Frequency');
title('Pilot Retransmission Frequency');
legend show;

% Helper function for AR simulation
function [retrans_freq, ber] = simulateAR(channel_values, use_coding, block_size, ...
    var_order, mse_threshold, num_initial_pilots, EbNo_dB, ldpc)
    
    % Initialize variables
    current_idx = 1;
    total_bits = 0;
    wrong_bits = 0;
    pilot_bits = 0;
    
    % Calculate maximum possible advance in one iteration
    initial_bits = num_initial_pilots * block_size;
    if use_coding
        max_block_size = 2*block_size;
    else
        max_block_size = block_size;
    end
    max_advance = max(max_block_size, initial_bits);
    
    % Get initial channel estimates
    channel_estimates = getInitialEstimates(channel_values, current_idx, ...
        num_initial_pilots, block_size, EbNo_dB);
    current_idx = current_idx + initial_bits;
    pilot_bits = pilot_bits + initial_bits;
    
    % Convert Eb/No to SNR
    if use_coding
        snr_db = EbNo_dB + 10*log10(1/2);  % Add coding rate adjustment
    else
        snr_db = EbNo_dB;
    end
    SNR = 10^(snr_db/10);
    noise_var = 1/(2*SNR);
    
    while current_idx + max_advance <= length(channel_values)
        % Train VAR model
        channel_data = [real(channel_estimates), imag(channel_estimates)];
        mdl = varm(2, var_order);
        mdl = estimate(mdl, channel_data);
        
        % Predict next channel value
        [Y_pred, ~] = forecast(mdl, 1, channel_data);
        htilde = complex(Y_pred(1), Y_pred(2));
        
        if use_coding
            [wrong_bits_block, hnought] = processCodedBlock(channel_values, current_idx, ...
                htilde, block_size, snr_db, ldpc);
            current_block_size = 2*block_size;
        else
            [wrong_bits_block, hnought] = processUncodedBlock(channel_values, current_idx, ...
                htilde, block_size, EbNo_dB);
            current_block_size = block_size;
        end
        
        % Check MSE
        mse = abs(htilde - hnought)^2;
        
        if mse > mse_threshold
            % Retransmit pilots
            channel_estimates = getInitialEstimates(channel_values, current_idx, ...
                num_initial_pilots, block_size, EbNo_dB);
            current_idx = current_idx + initial_bits;
            pilot_bits = pilot_bits + initial_bits;
        else
            % Update and continue
            channel_estimates = [channel_estimates(2:end); hnought];
            wrong_bits = wrong_bits + wrong_bits_block;
            total_bits = total_bits + current_block_size;
            current_idx = current_idx + current_block_size;
        end
    end
    
    % Calculate retransmission frequency (should be ≤ 1)
    retrans_freq = pilot_bits / (pilot_bits + total_bits);
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
    
    % Generate random information bits - use ldpc.K for correct size
    info_bits = randi([0 1], ldpc.K, 1);  % Changed from block_size to ldpc.K
    
    % Encode using LDPC
    coded_bits = ldpc.encode_bits(info_bits);  % Will produce ldpc.N bits
    symbols = 1 - 2*coded_bits;
    
    % Get actual channel values for full coded block
    actual_channel = channel_values(current_idx:current_idx+ldpc.N-1);  % Changed to ldpc.N
    
    % Add noise
    noise = sqrt(noise_var/2) * (randn(ldpc.N,1) + 1j*randn(ldpc.N,1));  % Changed size to ldpc.N
    received = actual_channel.*symbols + noise;
    
    % Compute LLRs using predicted channel
    llr = 2*real(conj(htilde).*received)./(abs(htilde)^2 * noise_var);
    
    % Decode using LDPC
    decoded_bits = ldpc.decode_llr(llr, 50, true);
    
    % Count errors (only in information bits)
    wrong_bits = sum(decoded_bits(1:ldpc.K) ~= info_bits);
    
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

