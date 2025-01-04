format long infsup
close all

INTERVAL_MODE = 1;

N = 1E+6; % Number of subintervals for s; 1e+5: NG, 1e+6:ok
s = linspace(0, 1, N+2);
s_intervals = I_infsup(s(1:end-2), s(3:end));

% X_list = [infsup(0,1),infsup(1.6,2.2),infsup(2.6,3.2),infsup(3.5,4)];
% X_list = [infsup(2.6,3.1),infsup(3.5,4)];
X_list = [infsup(2.6,3.1)]; % N=1e+5
X_list = [infsup(3.5,4)]; % N=1e+6

for j  = 1:length(X_list)
    X = X_list(j); % Search for non-zero kappa interval in this range
    
    % Resume processing from the next index
    for i = 1:N
        tic
        disp(i)
        s_val = s_intervals(i);
    
        if s_val < 0.9
            f = @(x)(func_left_hand_side(s_val, x));
        else
            sigma = (1-s_val)^(1/intval(3));
            f = @(x)(func_left_hand_side_singularity(s_val, x));
        end
    
        % Find the first three zeros of the function within the interval X
        no_sign_changes = find_non_zeros(f, X, 100);
        
        if any(no_sign_changes == 0)
            error('Sign changes detected in f_values.');
        end

        estimated_time = (N-i)*toc/60
    end
end

disp 'no sign changes for the given intervals'
