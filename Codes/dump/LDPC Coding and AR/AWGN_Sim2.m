% MATLAB Script: BER vs Eb/No for AWGN channel with LDPC
clear all;

% Parameters
block_length = 204;
info_length = 102;
rate = 1/2;                  % Code rate
max_iter = 50;              % Maximum number of decoding iterations
EbNo_dB = -1:0.2:5.4;      % Eb/No values in dB
num_frames = 100;           % Number of frames per Eb/No point
n_0 = 1/2;

%H = readAlist("code1.txt");
%H = sparse(logical(H));
rate = 1/2;
block_length = 1296;
info_length = block_length * rate;
[cfgLDPCEnc, cfgLDPCDec] = generateConfigLDPC(rate,block_length, 'decoderAlgo','norm-min-sum');

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
     
        % Encode using custom LDPC code
        coded_bits = ldpcEncode(info_bits, cfgLDPCEnc);
       
        % BPSK modulation (1 -> +1, 0 -> -1)
        symbols = 1 - 2 * coded_bits;
        
        % Generate and add AWGN noise
        noise = sqrt(noise_var) * randn(block_length, 1);
        received = symbols + noise;
        
        % Compute LLR
        llr = 2 * received/noise_var;
        
        % Decode using LDPC
        decoded_bits = ldpcDecode(llr,cfgLDPCDec,max_iter);
        
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

function [H] = readAlist(fname) 
% reads binary parity check matrix in "alist" format from file FNAME and 
% converts it to sparse matrix used in MATLAB routines. 
% This is an interface to matrices at http://wol.ra.phy.cam.ac.uk/mackay/codes/ 
% 
% Example 
%        [H] = alist2sparse('A');   % A is the ascii file in alist format 
 
 
%   Copyright (c) 1999 by Igor Kozintsev igor@ifp.uiuc.edu 
%   $Revision: 1.1 $  $Date: 2000/03/23 $ Bug fixed by Hatim Behairy 
 
fid = fopen(fname); 
n = fscanf(fid,'%d',1); 
m = fscanf(fid,'%d',1); 
maxinrow = fscanf(fid,'%d',1);  
junk = fscanf(fid,'%d',1); % no need 
num = fscanf(fid,'%d',[1 n]); % number of elements in rows 
 
num2(1:n)=maxinrow;     
junk = fscanf(fid,'%d',[1 m]); % no need 
 
position = zeros(n,maxinrow); 
for i=1:n 
   for j=1:num2(i)     
      position(i,j) = fscanf(fid,'%d',1); 
   end 
end 
 
ii = zeros(1,sum(num)); 
jj = ii; 
k = 1; 
for i=1:n 
      for j=1:num(i) 
      jj(k) = i; 
      ii(k) = position(i,j); 
      ss = 1; 
      k = k+1 ;  
   end 
end 
H = sparse(ii,jj,ss,m,n); 
fclose(fid);
end