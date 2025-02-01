% MATLAB Script: BER vs Eb/No for AWGN channel with LDPC
clear all;

% Parameters
block_length = 96;
info_length = 48;                 % Code rate
max_iter = 100;              % Maximum number of decoding iterations
EbNo_dB = -1:0.2:5.4;      % Eb/No values in dB
num_frames = 100;           % Number of frames per Eb/No point
n_0 = 1/2;                  % Noise power spectral density

H = N96M48();
[M, N, K, vn_degree, cn_degree, P, H_row_one_absolute_index, H_comlumn_one_relative_index, vn_distribution, cn_distribution] = H_matrix_process(H);
%cfgLDPCEnc = ldpcEncoderConfig;
%decodercfg = ldpcDecoderConfig(cfgLDPCEnc);
rate = K/N;
max_runs = 1e8;
resolution = 1e5;

% Convert Eb/No to SNR
snr_db_vec = EbNo_dB + 10*log10(5/6); 

% Preallocate BER arrays
BER_coded = zeros(size(EbNo_dB));
BER_uncoded = zeros(size(EbNo_dB));

% Loop over Eb/No values
for i = 1:length(EbNo_dB)
    fprintf('Simulating Eb/No = %.1f dB\n', EbNo_dB(i));
    SNR_linear = 10^(-snr_db_vec(i)/20);
    noise_var = 1 / (sqrt(2 * rate)*SNR_linear);
  
    bit_errors_coded = 0;
    total_bits_coded = 0;
    bit_errors_uncoded = 0;
    total_bits_uncoded = 0;
    
    for frame = 1:num_frames
        %% Coded Transmission
        % Generate random information bits
        info_bits = randi([0 1],cfgLDPCEnc.NumInformationBits,1);
        parity_check_bits = mod(P * info_bits, 2);
        coded_bits = [info_bits; parity_check_bits];
        
        % Encode the bits
        % coded_bits = ldpcEncode(info_bits_coded, cfgLDPCEnc);
       
        % BPSK modulation (1 -> +1, 0 -> -1)
        symbols = 1 - 2 * coded_bits;
        
        % Generate and add AWGN noise
        noise = noise_var * randn(block_length, 1);
        received = symbols + noise;
        
        % Compute LLR
        llr = 2*received/noise_var^2;
        
        % Decode using LDPC
        %[decoded_bits, actualnumiter, finalchecks] = ldpcDecode(llr, decodercfg, max_iter);
        [decoded_bits, ~] = bp_decoder(llr, H, max_iter);
        
        % Count bit errors for coded transmission
        bit_errors_coded = bit_errors_coded + sum(decoded_bits ~= info_bits);
        total_bits_coded = total_bits_coded + info_length;
        
        %% Uncoded Transmission
        SNR_linear = 10^(-EbNo_dB(i)/20);
        noise_var = 1 / (sqrt(2) * SNR_linear);

        % Generate random information bits
        info_bits_uncoded = rand(block_length, 1) < 0.5;
        
        % BPSK modulation
        symbols_uncoded = 1 - 2 * info_bits_uncoded;
        
        % Add AWGN noise
        noise_uncoded = noise_var * randn(block_length, 1);
        received_uncoded = symbols_uncoded + noise_uncoded;
        
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

function H_sparse = N96M48()
H = [34	93	26	50	8	49
32	76	42	72	6	67
30	54	14	80	38	81
46	88	8	90	47	55
13	84	38	69	26	59
41	92	30	88	13	76
44	68	29	55	21	89
35	70	11	74	34	88
1	55	20	94	37	63
40	80	32	77	29	52
4	52	2	83	1	50
7	81	27	59	46	64
5	69	4	92	9	93
45	74	6	52	14	92
38	67	15	71	27	74
21	65	31	57	7	73
3	62	37	54	33	71
23	94	45	84	4	95
12	89	23	65	32	58
15	77	35	58	42	54
17	60	21	75	5	87
22	64	5	78	20	79
42	83	33	79	12	85
16	50	13	62	40	86
2	63	44	70	28	57
26	91	46	56	23	82
33	57	25	60	15	80
27	72	9	61	48	75
48	71	47	81	16	77
9	82	17	68	45	62
11	51	43	73	22	72
31	90	28	82	3	96
19	73	10	63	30	56
43	59	3	67	18	70
37	85	41	86	17	66
39	86	18	95	10	60
10	58	16	91	31	51
18	56	40	76	44	78
6	66	7	49	11	53
29	78	36	93	19	65
20	96	24	89	36	83
25	49	39	64	2	90
36	61	22	53	35	91
28	95	34	85	25	61
47	53	19	87	41	84
24	79	48	51	39	69
8	75	1	96	43	68
14	87	12	66	24	94];
H_sparse = zeros(size(H, 1), max(max(H)));
for m = 1 : size(H, 1)
    ones_index = H(m, :);
    ones_index = ones_index(ones_index ~= 0);
    H_sparse(m, ones_index) = 1;
end
end

% LDPC Matrix Processing Functions
% This script contains utility functions for processing LDPC (Low-Density Parity-Check) matrices

function H_output = Gaussian_Elimination(H_input)
    % Performs Gaussian elimination on a binary matrix
    % Input: H_input - binary matrix
    % Output: H_output - row-reduced matrix
    
    N = size(H_input, 2);
    M = size(H_input, 1);
    level = 0;
    
    for i = N : -1 : 1  % Gaussian elimination
        next_column = 0;
        for j = 1 : M - level
            if H_input(j, i) == 1
                tmp = H_input(j, :);
                H_input(j, :) = H_input(M - level, :);
                H_input(M - level, :) = tmp;
                level = level + 1;
                break;
            else
                if j == M - level
                    next_column = 1;
                end
            end
        end
        
        if next_column == 1
            continue;
        end
        
        % Clear 1's upward
        for k = 1 : M - level
            if H_input(k, i) == 1
                H_input(k, :) = mod(H_input(k, :) + H_input(M - level + 1, :), 2);
            end
        end
        
        % Clear 1's downward
        for k = M - level + 2 : M
            if H_input(k, i) == 1
                H_input(k, :) = mod(H_input(k, :) + H_input(M - level + 1, :), 2);
            end
        end
        
        if level == M
            break;
        end
    end
    
    % Remove all-zero rows
    rank_loss = 0;
    while(all(H_input(1, :) == 0))
        H_input(1, :) = [];
        rank_loss = rank_loss + 1;
    end
    
    % Display rank information
    if rank_loss ~= 0
        disp(['There is redundancy in H. The number of dependent rows is ' num2str(rank_loss) '.'])
    else
        disp('The initial H is full row rank.')
    end
    
    H_output = H_input;
end

function [H, P] = H_systematic(H_simplified_stair, H_original)
    % Convert H matrix to systematic form
    % Input: 
    %   H_simplified_stair - Simplified staircase matrix
    %   H_original - Original parity-check matrix
    % Output:
    %   H - Systematic form of the matrix
    %   P - Submatrix used for parity check bits
    
    N = size(H_simplified_stair, 2);
    K = size(H_simplified_stair, 1);  % Number of parity check bits
    
    H = zeros(size(H_original, 1), N, 'uint8');
    I_location = zeros(K, 1);
    
    % Find identity column locations
    for i = N : -1 : 1
        index = find(H_simplified_stair(:, i));
        if (length(index) == 1) && (I_location(index) == 0)
           I_location(index) = i;
        end
    end
    
    % Permute columns if necessary
    if any(I_location ~= (N - K + 1 : N)')
        disp('Note: Permuting columns of H_original to get systematic LDPC codes.')
        P_location = set_minus(1 : N, I_location);
        H(:, 1 : N - K) = H_original(:, P_location);
        H(:, N - K + 1 : N) = H_original(:, I_location);
        P = H_simplified_stair(:, P_location);
    else
        disp('No Column Permutations.')
        H = H_original;
        P = H_simplified_stair(:, 1 : N - K);
    end
end

function [H_row_one_absolute_index, H_column_one_relative_index, vn_degree, cn_degree, vn_distribution, cn_distribution] = extract_H_structure(H)
    % Extract structural information from H matrix
    % Input: H - Parity-check matrix
    % Output: Various structural parameters of the matrix
    
    % Variable node (VN) and check node (CN) degrees
    vn_degree = sum(H);
    cn_degree = sum(H, 2);
    
    M = size(H, 1);
    N = size(H, 2);
    
    % Find absolute indices of 1's in each row
    H_row_one_absolute_index = zeros(M, max(cn_degree));
    for i = 1 : M
        cnt = 1;
        for j = 1 : N
            if H(i, j) == 1
                H_row_one_absolute_index(i, cnt) = j;
                cnt = cnt + 1;
            end
        end
    end
    
    % Find relative indices of 1's in each column
    H_column_one_relative_index = zeros(M, max(cn_degree));
    for i = 1 : M
        for j = 1 : cn_degree(i)
            H_column_one_relative_index(i, j) = 1 + sum(H(1 : i - 1, H_row_one_absolute_index(i, j)));
        end
    end
    
    % Variable node degree distribution
    vn_distinct_degree_tmp = zeros(1, max(vn_degree));
    for i = 1 : N
        vn_distinct_degree_tmp(vn_degree(i)) = vn_distinct_degree_tmp(vn_degree(i)) + 1;
    end
    distinct_degree = find(vn_distinct_degree_tmp);
    degree_distribution = vn_distinct_degree_tmp(distinct_degree)/N;
    vn_distribution = num2str([distinct_degree; degree_distribution]);
    
    % Check node degree distribution
    cn_distinct_degree_tmp = zeros(1, max(cn_degree));
    for i = 1 : M
        cn_distinct_degree_tmp(cn_degree(i)) = cn_distinct_degree_tmp(cn_degree(i)) + 1;
    end
    distinct_degree = find(cn_distinct_degree_tmp);
    degree_distribution = cn_distinct_degree_tmp(distinct_degree)/M;
    cn_distribution = num2str([distinct_degree; degree_distribution]);
end

function C = set_minus(A, B)
    % Set difference operation
    % Returns elements in A that are not in B
    
    for i = 1 : length(A)
        if ismember(A(i), B)
            A(i) = -1;
        end
    end
    C = A(A > 0);
end

function [M, N, K, vn_degree, cn_degree, P, H_row_one_absolute_index, H_column_one_relative_index, vn_distribution, cn_distribution] = H_matrix_process(H_input)
    % Process LDPC parity-check matrix
    % Input: H_input - Original sparse binary matrix
    % Output: Various matrix parameters and structural information
    
    % Matrix dimensions
    M = size(H_input, 1);
    N = size(H_input, 2);
    
    % Gaussian elimination
    H_simplified_stair = Gaussian_Elimination(H_input);
    
    % Number of information bits
    K = N - size(H_simplified_stair, 1);
    
    % Convert to systematic form
    [H_column_permuted, P] = H_systematic(H_simplified_stair, H_input);
    
    % Display and convert P to double for calculations
    disp(['The data-class of P is ' class(P)])
    P = double(P);
    
    % Extract matrix structure
    [H_row_one_absolute_index, H_column_one_relative_index, vn_degree, cn_degree, vn_distribution, cn_distribution] = extract_H_structure(H_column_permuted);
end