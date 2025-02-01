clear all;
close all;

% Parameters
block_length = 1296; % Block length
rate = 1/2; % Code rate
max_iter = 50; % Maximum number of decoding iterations
EbNo_dB = -1:0.5:20; % Eb/No values in dB
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
%Convert Eb/No to SNR
snr_db_vec = EbNo_dB + 10*log10(rate); 

% Track current position in channel values
current_index = 1;

% Loop over Eb/No values
for i = 1:length(EbNo_dB)
    
    fprintf('Simulating Eb/No = %.1f dB\n', EbNo_dB(i));
    bit_errors_uncoded2 = 0;
    total_bits_uncoded2 = 0;
    frame = 0;
    %{
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
    %}
    while (frame < max_frames && current_index + 2*block_length <= length(channel_values))
        frame = frame + 1;
        
        %{
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
        noise_coded = sqrt(noise_var) * (randn(ldpc.N, 1) + 1j * randn(ldpc.N, 1));
        
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
        noise_uncoded = sqrt(noise_var) * (randn(block_length, 1) + 1j * randn(block_length, 1));
        
        % Pass through channel with noise
        received_uncoded = channel_uncoded .* symbols_uncoded + noise_uncoded/sqrt(SNR_linear);
        
        % Hard decision decoding
        decoded_uncoded = real(received_uncoded ./ channel_uncoded) < 0;
        
        % Count bit errors for uncoded transmission
        bit_errors_uncoded1 = bit_errors_uncoded1 + sum(decoded_uncoded ~= info_bits_uncoded);
        total_bits_uncoded1 = total_bits_uncoded1 + block_length;
        %}
        %%Uncoded Transmission 2
        % Generate random information bits
        SNR_linear = 10^(EbNo_dB(i)/10);
        noise_var = 1 / (2 * SNR_linear);
        info_bits_uncoded = randi([0, 1], block_length, 1);

        % BPSK modulation
        symbols_uncoded = 1 - 2 * info_bits_uncoded;

        % Generate Rayleigh channel and complex noise
        channel_uncoded = sqrt(0.5) * (randn(block_length, 1) + 1j * randn(block_length, 1));
        noise_uncoded = sqrt(noise_var) * (randn(block_length, 1) + 1j * randn(block_length, 1));

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
                frame, bit_errors_coded/total_bits_coded, bit_errors_uncoded/total_bits_uncoded);
        end
    end
    
    % Reset channel index for next SNR point
    current_index = 1;
    
    % Calculate BER for the current Eb/No value
    %{
    BER_coded(i) = bit_errors_coded / total_bits_coded;
    BER_uncoded1(i) = bit_errors_uncoded1 / total_bits_uncoded1;
    %}
    BER_uncoded2(i) = bit_errors_uncoded2 / total_bits_uncoded2;
end

%}

% Now initialize arrays for AR simulation
BER_ar_uncoded = zeros(size(EbNo_dB));
BER_ar_coded = zeros(size(EbNo_dB));
retrans_freq_uncoded = zeros(size(EbNo_dB));
retrans_freq_coded = zeros(size(EbNo_dB));


% AR parameters
var_order = 25;
num_initial_pilots = 200;
mse_threshold = 0.1;

%Convert Eb/No to SNR
snr_db_vec = EbNo_dB + 10*log10(rate);

% Run AR simulations
for i = 1:length(EbNo_dB)
    fprintf('Simulating AR for Eb/No = %.1f dB\n', EbNo_dB(i));
    
    [retrans_freq_uncoded(i), BER_ar_uncoded(i)] = simulateARuncoded(channel_values, var_order, mse_threshold, num_initial_pilots, EbNo_dB);
    
    [retrans_freq_coded(i), BER_ar_coded(i)] = simulateARcoded(channel_values, var_order, mse_threshold, num_initial_pilots, snr_db_vec);
end


% Plot all results
figure(1);
%semilogy(EbNo_dB, BER_coded, '-o', 'LineWidth', 2, 'DisplayName', 'Coded Jakes Channel');
hold on;
%semilogy(EbNo_dB, BER_uncoded1, '-s', 'LineWidth', 2, 'DisplayName', 'Jakes Channel');
semilogy(EbNo_dB, BER_uncoded2, '-s', 'LineWidth', 2, 'DisplayName', 'Rayleigh Channel');
semilogy(EbNo_dB, BER_ar_uncoded, '-v', 'LineWidth', 2, 'DisplayName', 'AR with Jakes Channel - Coded');
semilogy(EbNo_dB, BER_ar_coded, '-d', 'LineWidth', 2, 'DisplayName', 'AR with Jakes Channel - Uncoded');
grid on;
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs Eb/No for Coded and Uncoded Transmission');
legend show;

figure(2);
plot(EbNo_dB, retrans_freq_uncoded, '-o', 'DisplayName', 'AR with Jakes Channel - Uncoded');
hold on;
plot(EbNo_dB, retrans_freq_coded, '-s', 'DisplayName', 'AR with Jakes Channel - Coded');
grid on;
xlabel('Eb/No (dB)');
ylabel('Retransmission Frequency');
title('Pilot Retransmission Frequency');
legend show;


function [retrans_freq, ber] = simulateARuncoded(channel_values, var_order, mse_threshold, num_initial_pilots, EbNo_dB)
    % Initialize variables
    current_idx = 1;
    total_bits = 0;
    wrong_bits = 0;
    pilot_bits = 0;
    pilot_block_size = 1;
    block_size = 1;

    % Calculate maximum possible advance in one iteration
    initial_bits = num_initial_pilots * block_size;
    max_block_size = 200;
    max_advance = 200;
    
    % Get initial channel estimates
    initial_channel_estimates = getInitialEstimates(channel_values, current_idx, num_initial_pilots, pilot_block_size, EbNo_dB);
    % Simulation variables
    current_channel_estimates = initial_channel_estimates;
    current_idx = current_idx + initial_bits;
    pilot_bits = pilot_bits + initial_bits;
    
    snr_db = EbNo_dB;
    SNR = 10.^(snr_db/10);
    noise_var = 1/(2*SNR);
    
    while current_idx + 2*max_advance <= length(channel_values)
        % Train VAR model
        channel_data = [real(current_channel_estimates), imag(current_channel_estimates)];
        mdl = varm(2, var_order);
        mdl = estimate(mdl, channel_data);
        
        % Predict next channel value
        [Y_pred, ~] = forecast(mdl, 10, channel_data);
        htilde = complex(Y_pred(:, 1), Y_pred(:, 2));


        info_bits = randi([0, 1], 10, 1);
        symbols = 1 - 2 * info_bits;
        % Get actual channel values
        actual_channels = channel_values(current_idx:current_idx+9);
        current_idx = current_idx + 10;
        % Add noise
        SNR = 10^(EbNo_dB/10);
        noise_var = 1 / (2 * SNR);
        noise = sqrt(noise_var) * (randn(10, 1) + 1j * randn(10, 1));
        received = actual_channels .* symbols + noise/sqrt(SNR);

        predicted_bits = real(received./htilde) < 0;

        hnought = received./predicted_bits;

        wrong_bits = wrong_bits + sum(predicted_bits ~= info_bits);
        total_bits = total_bits + 10;
     
        % Check MSE
        mse = abs(htilde - hnought)^2;
        if mse > mse_threshold
            % Retransmit 200 pilots
            initial_channel_estimates = getInitialChannelEstimates(channel_values, current_idx, num_initial_pilots, block_length, EbNo_dB);
            current_idx = current_idx + num_initial_pilots * block_length;
            pilot_bits = pilot_bits + num_initial_pilots * block_length;
            current_channel_estimates = initial_channel_estimates;
        else
            current_channel_estimates = [current_channel_estimates(11:end); hnought];            
        end
    end
    
    % Calculate retransmission frequency and BER
    retrans_freq = pilot_bits / (pilot_bits + total_bits);
    ber = wrong_bits / total_bits;
end


function [retrans_freq, ber] = simulateARcoded(channel_values, var_order, mse_threshold, num_initial_pilots, EbNo_dB)
    % Initialize variables
    current_idx = 1;
    total_bits = 0;
    wrong_bits = 0;
    pilot_bits = 0;
    pilot_block_size = 1;
    block_size = 1;

    snr_db = EbNo_dB + 10*log10(rate);
    SNR = 10.^(snr_db/10);
    noise_var = 1/(2*SNR);

    % Calculate maximum possible advance in one iteration
    initial_bits = num_initial_pilots * block_size;
    max_block_size = 200;
    max_advance = 200;
    
    % Get initial channel estimates
    initial_channel_estimates = getInitialEstimates(channel_values, current_idx, num_initial_pilots, pilot_block_size, snr_db);
    % Simulation variables
    current_channel_estimates = initial_channel_estimates;
    current_idx = current_idx + initial_bits;
    pilot_bits = pilot_bits + initial_bits;
    
    while current_idx + 2*max_advance <= length(channel_values)
        % Train VAR model
        channel_data = [real(current_channel_estimates), imag(current_channel_estimates)];
        mdl = varm(2, var_order);
        mdl = estimate(mdl, channel_data);
        
        % Predict next channel value
        [Y_pred, ~] = forecast(mdl, 1296, channel_data);
        htilde = complex(Y_pred(:, 1), Y_pred(:, 2));


        info_bits = randi([0, 1], 648, 1);
        % Encode the bits
        codeword = ldpc.encode_bits(info_bits_coded);
        symbols_coded = 1 - 2 * codeword;
        % Get channel values for this frame
        channel_coded = channel_values(current_index:current_index + 1295);
        current_index = current_index + 1296;       
        % Generate complex noise
        noise_coded = sqrt(noise_var/2) * (randn(ldpc.N, 1) + 1j * randn(ldpc.N, 1));
        % Pass through channel with noise
        received_coded = channel_coded .* symbols_coded + noise_coded/sqrt(snr);
        % Compute LLRs
        llr_coded = real(2 * conj(htilde) .* received_coded ./ (abs(htilde).^2 * noise_var));
        % Decode using LDPC
        decoded_codeword = ldpc.decode_llr(llr_coded, max_iter, true);

        predicted_bits = ldpc.encode_bits(decoded_codeword);

        hnought = received_coded./predicted_bits;

        wrong_bits = wrong_bits + sum(decoded_codeword ~= info_bits);
        total_bits = total_bits + 648;
     
        % Check MSE
        mse = abs(htilde - hnought)^2;
        if mse > mse_threshold
            % Retransmit 200 pilots
            initial_channel_estimates = getInitialChannelEstimates(channel_values, current_idx, num_initial_pilots, block_length, EbNo_dB);
            current_idx = current_idx + num_initial_pilots * block_length;
            pilot_bits = pilot_bits + num_initial_pilots * block_length;
            current_channel_estimates = initial_channel_estimates;
        else
            current_channel_estimates = [hnought(end-200:end)];            
        end
    end
    
    % Calculate retransmission frequency and BER
    retrans_freq = pilot_bits / (pilot_bits + total_bits);
    ber = wrong_bits / total_bits;
end

function channel_estimates = getInitialEstimates(channel_values, current_idx, num_pilots, block_size, EbNo_dB)
    % Initialize channel estimates array
    channel_estimates = zeros(num_pilots, 1);
    
    % Convert Eb/No to SNR
    SNR = 10.^(EbNo_dB/10);
    
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



