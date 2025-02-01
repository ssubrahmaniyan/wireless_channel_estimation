% MATLAB Script: BER vs Eb/No for AWGN channel with BPSK using LDPC
clear all;

% Parameters
block_length = 1296;          % Block length
rate = 1/2;                  % Code rate
max_iter = 50;              % Maximum number of decoding iterations
EbNo_dB = -1:1:10;         % Eb/No values in dB
num_frames = 100;           % Number of frames per Eb/No point
n_0 = 1/2;                  % Noise power spectral density

% Initialize LDPC code
ldpc = LDPCCode(0, 0);
ldpc.load_wifi_ldpc(block_length, rate);
info_length = ldpc.K;

% Convert Eb/No to SNR
snr_db_vec = EbNo_dB + 10*log10(info_length/block_length); 

% Preallocate BER arrays
BER_coded = zeros(size(EbNo_dB));
BER_uncoded = zeros(size(EbNo_dB));

% Loop over Eb/No values
for i = 1:length(EbNo_dB)
    fprintf('Simulating Eb/No = %.1f dB\n', EbNo_dB(i));
    SNR_linear = 10^(snr_db_vec(i)/10);
    noise_var = 1 / (2 * rate * SNR_linear);
    
    snr = 10^(snr_db_vec(i)/10);
    bit_errors_coded = 0;
    total_bits_coded = 0;
    bit_errors_uncoded = 0;
    total_bits_uncoded = 0;
    
    for frame = 1:num_frames
        %% Coded Transmission
        % Generate random information bits
        info_bits = rand(info_length, 1) < 0.5;
        
        % Encode the bits
        coded_bits = ldpc.encode_bits(info_bits);
       
        % BPSK modulation (1 -> +1, 0 -> -1)
        symbols = 1 - 2 * coded_bits;
        
        % Generate and add AWGN noise
        noise = sqrt(1/2) * randn(block_length, 1);
        received = symbols + noise/sqrt(snr);
        
        % Compute LLR
        llr = 2 * received/noise_var;
        
        % Decode using LDPC
        [decoded_bits, ~] = ldpc.decode_llr(llr, max_iter, true);
        
        % Count bit errors for coded transmission
        bit_errors_coded = bit_errors_coded + sum(decoded_bits(1:info_length) ~= info_bits);
        total_bits_coded = total_bits_coded + info_length;
        
        %% Uncoded Transmission
        SNR_linear = 10^(EbNo_dB(i)/10);
        noise_var = 1 / (2 * SNR_linear);

        % Generate random information bits
        info_bits_uncoded = rand(block_length, 1) < 0.5;
        
        % BPSK modulation
        symbols_uncoded = 1 - 2 * info_bits_uncoded;
        
        % Add AWGN noise
        noise_uncoded = sqrt(n_0) * randn(block_length, 1);
        received_uncoded = symbols_uncoded + noise_uncoded/sqrt(SNR_linear);
        
        % Hard decision decoding
        decoded_uncoded = received_uncoded < 0;
        
        % Count bit errors for uncoded transmission
        bit_errors_uncoded = bit_errors_uncoded + sum(decoded_uncoded ~= info_bits_uncoded);
        total_bits_uncoded = total_bits_uncoded + block_length;
    end
    
    % Calculate BER for the current Eb/No value
    BER_coded(i) = bit_errors_coded / total_bits_coded;
    BER_uncoded(i) = bit_errors_uncoded / total_bits_uncoded;
end

% Plot BER vs Eb/No
figure;
semilogy(EbNo_dB, BER_coded, '-o', 'LineWidth', 2, 'DisplayName', 'LDPC Coded BPSK');
hold on;
semilogy(EbNo_dB, BER_uncoded, '-s', 'LineWidth', 2, 'DisplayName', 'Uncoded BPSK');
grid on;
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs Eb/No for Coded and Uncoded BPSK over AWGN Channel');
legend show;