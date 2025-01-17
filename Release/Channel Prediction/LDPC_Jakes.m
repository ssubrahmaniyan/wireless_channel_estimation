% Parameters
block_length = 648; % Block length
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
                frame, bit_errors_coded/total_bits_coded, bit_errors_uncoded/total_bits_uncoded);
        end
    end
    
    % Reset channel index for next SNR point
    current_index = 1;
    
    % Calculate BER for the current Eb/No value
    BER_coded(i) = bit_errors_coded / total_bits_coded;
    BER_uncoded1(i) = bit_errors_uncoded1 / total_bits_uncoded1;
    BER_uncoded2(i) = bit_errors_uncoded2 / total_bits_uncoded2;
end

% Plot BER vs Eb/No
figure;
semilogy(EbNo_dB, BER_coded, '-o', 'LineWidth', 2, 'DisplayName', 'Coded');
hold on;
semilogy(EbNo_dB, BER_uncoded1, '-s', 'LineWidth', 2, 'DisplayName', 'Uncoded 1');
semilogy(EbNo_dB, BER_uncoded2, '-s', 'LineWidth', 2, 'DisplayName', 'Uncoded 2');
grid on;
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs Eb/No for Coded and Uncoded Transmission (Jakes, BPSK)');
legend show;