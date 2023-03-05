% gauss_elim_P2p1
function [x] = gauss_elim_P2p1(A, b)
format compact;

% =========================================================================
% ================================ Main ===================================
% =========================================================================

% Is matrix A square?
if is_square(A) % If true, compute x.
    aug = [A b];  % Augment A and b.
    aug = row_echelon(aug); % Get equivalent GE coefficent matrix in row-echelon form.
    x = back_substitution(aug); % Find x using back substitution.
else % Else, error message.
    msg = 'Matrix A is not a square matrix!';
    error(msg);
end

% =========================================================================
% ======================== Helper Sub Functions ===========================
% =========================================================================

% Checking if A is a square matrix.
    function boolean = is_square(A)
        s = size(A); % Get number of rows and columns.
        if s(1) ~= s(2) % Are num of rows equal to num of cols?
            boolean = false; % No, retrun false.
        else
            boolean = true; % Yes, return true.
        end
    end

% Row echelon procedure
    function aug = row_echelon(aug)
        [row, col] = size(aug); % Get number of rows and columns.
        N = col-1; % Num of columns for A.
        aug = near_zero_pivot_check(1 , 1, aug); % Check magnitude of first pivot.
        for k = 1:N
            for i = k+1:row
                aug = near_zero_pivot_check(i , k, aug); % Check magnitude of ith kith pivot.
                lambda = aug(i, k)/aug(k, k);
                aug(i, k:N) = aug(i, k:N) - lambda * aug(k, k:N);
                aug(i,col) = aug(i,col) - lambda * aug(k,col);
            end
        end
    end

% Check if pivot is near zero.
    function aug = near_zero_pivot_check(i , k, aug)
        near_zero = 0.02;
        abs_pivot = abs(aug(i,k));
        if abs_pivot < near_zero % Is pivot less than 0.02?
            largest_pivot_index = search_for_largest_pivot(i, k, aug); % Find row index of largest pivot.
            aug = swap(i, largest_pivot_index, aug); % Swap rows at index i and largest_pivot_index.
        end
    end

% Finds row index with largest pivot
    function largest_pivot_index = search_for_largest_pivot(i, k, aug)
        r = length(aug) - 1; % Row num = (num of cols - 1).
        largest_pivot_index = i; % Assume i is the largest pivot.
        for row = i+1:r
            if aug(largest_pivot_index,k) < aug(row,k) % Is largest pivot less than current pivot?
                largest_pivot_index = row; % Set largest pivot index to current row index.
            end
        end
    end

% Swaps rows
    function aug =  swap(i, largest_index, aug)
        temp = aug(i,:);
        aug(i,:) = aug(largest_index,:); 
        aug(largest_index,:) = temp; 
    end

% Back Substitution
    function x = back_substitution(aug)
        [r, N] = size(aug);
        A_ = aug(:, 1:N-1);
        b_ = aug(:,N);
        x = b_;
        x(r) = x(r)/A_(r,r);
        for i = r:-1:1
            sum = 0;
            for j = i+1:(N-1)
                sum = sum + A_(i,j)*x(j);
            end
            x(i) = (b_(i) - sum) / A_(i,i);
        end
    end
end