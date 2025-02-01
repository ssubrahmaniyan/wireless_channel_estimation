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