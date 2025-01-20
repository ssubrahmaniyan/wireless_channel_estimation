% MATLAB Script: BER vs Eb/No for AWGN channel with BPSK using LDPC and Rayleigh Channel

% Parameters
block_length = 1296;  % Block length
rate = 1/2;          % Code rate
max_iter = 100;       % Maximum number of decoding iterations
EbNo_dB = -1:1:30;  % Eb/No values in dB
num_frames = 200;    % Number of frames per Eb/No point

% Initialize LDPC code
ldpc = LDPCCode(block_length, block_length * rate);
ldpc.load_wifi_ldpc(block_length, rate);

% Preallocate BER arrays
BER_coded = zeros(size(EbNo_dB));
BER_uncoded = zeros(size(EbNo_dB));

% Convert Eb/No to SNR
snr_db_vec = EbNo_dB + 10*log10(info_length/block_length); 

% Loop over Eb/No values
for i = 1:length(EbNo_dB)
    fprintf('Simulating Eb/No = %.1f dB\n', EbNo_dB(i));

    bit_errors_coded = 0;
    total_bits_coded = 0;

    bit_errors_uncoded = 0;
    total_bits_uncoded = 0;

    for frame = 1:num_frames
        SNR_linear = 10^(snr_db_vec(i)/10);
        noise_var = 1 / (2 * rate * SNR_linear);
        snr = 10^(snr_db_vec(i)/10);
        %% Coded Transmission
        % Generate random information bits
        info_bits_coded = randi([0, 1], ldpc.K, 1);

        % Encode the bits
        codeword = ldpc.encode_bits(info_bits_coded);

        % BPSK modulation
        symbols_coded = 1 - 2 * codeword;

        % Generate Rayleigh channel and complex noise
        channel_coded = sqrt(0.5) * (randn(ldpc.N, 1) + 1j * randn(ldpc.N, 1));
        noise_coded = sqrt(1 / 2) * (randn(ldpc.N, 1) + 1j * randn(ldpc.N, 1));

        % Pass through channel with noise
        received_coded = channel_coded .* symbols_coded + noise_coded/sqrt(snr);

        % Compute LLRs
        llr_coded = real(2 * conj(channel_coded) .* received_coded ./ (abs(channel_coded).^2 * noise_var));

        % Decode using LDPC
        decoded_codeword = ldpc.decode_llr(llr_coded, max_iter, true);

        % Count bit errors for coded transmission
        bit_errors_coded = bit_errors_coded + sum(decoded_codeword(1:ldpc.K) ~= info_bits_coded);
        total_bits_coded = total_bits_coded + ldpc.K;

        %% Uncoded Transmission
        SNR_linear = 10^(EbNo_dB(i)/10);
        noise_var = 1 / (2 * SNR_linear);
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
        bit_errors_uncoded = bit_errors_uncoded + sum(decoded_uncoded ~= info_bits_uncoded);
        total_bits_uncoded = total_bits_uncoded + block_length;
    end

    % Calculate BER for the current Eb/No value
    BER_coded(i) = bit_errors_coded / total_bits_coded;
    BER_uncoded(i) = bit_errors_uncoded / total_bits_uncoded;
end

% Plot BER vs Eb/No
figure;
semilogy(EbNo_dB, BER_coded, '-o', 'LineWidth', 2, 'DisplayName', 'Coded');
hold on;
semilogy(EbNo_dB, BER_uncoded, '-s', 'LineWidth', 2, 'DisplayName', 'Uncoded');
grid on;
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs Eb/No for Coded and Uncoded Transmission (Rayleigh, BPSK)');
legend show;
