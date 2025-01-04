format long infsup
INTERVAL_MODE=1;

function val = compute_non_zeros(s_list)
    val =[];
    for i = 1:len(s_list)
        s_val = s_list(i);
    
        if s_val < 0.9
            f = @(x)(func_left_hand_side(s_val, x));
        else
            sigma = (1-s_val)^(1/intval(3));
            f = @(x)(func_left_hand_side_singularity(s_val, x));
        end
    
        f0 = f(inf(X));
    
        % Find the first three zeros of the function within the interval X
        non_zeros_for_each_s = find_non_zeros(f, f0, X, 200);
        lower_bounds_zeros = inf(non_zeros_for_each_s);
        [~, idx] = sort(lower_bounds_zeros);
        non_zeros_for_each_s = non_zeros_for_each_s(idx);
        lowerBound = min(inf(non_zeros_for_each_s));
        upperBound = max(sup(non_zeros_for_each_s));
        non_zeros_for_each_s = infsup(lowerBound, upperBound);
        
        % Extract lower and upper bounds of each interval
        zeros_for_each_s_L = I_inf(non_zeros_for_each_s);
        zeros_for_each_s_U = I_sup(non_zeros_for_each_s);
    
        % Interleave lower and upper bounds for clarity: [L1, U1, L2, U2, ...]
        M = length(zeros_for_each_s_L);
        intervals_flat = zeros(1, 2*M);
        intervals_flat(1:2:end) = zeros_for_each_s_L; % Place lower bounds at odd indices
        intervals_flat(2:2:end) = zeros_for_each_s_U; % Place upper bounds at even indices
    
        % Combine the index and the intervals into one row
        result_for_each_s = [i, intervals_flat]
    end

end

N1=6;
N2=15; % 15:ok, 14:NG, 13:NG
x0 = infsup(0,4);

% If the CSV file exists, read the last processed index
csv_file = 'zeros.csv';
if exist(csv_file, 'file')
    data = readmatrix(csv_file);
    if ~isempty(data)
        last_i = data(end, 1); % Retrieve the last processed index from the CSV file
    else
        last_i = 0; % Initialize to 0 if the file is empty
    end
else
    last_i = 0; % Initialize to 0 if the file does not exist
end

bi_N1 = 2^(-N1);
bi_N2 = 2^(-N2);

%treat the first some intervals differently because dF/Dk vanishses at s=0
i=1; i_=2^6;
s = hull(i-1,i)*bi_N1
f = @(x)(func_left_hand_side(s,x)); % Use singularity function if s < 0.05

options = verifynlssallset( ...
    'Boxes', 2^10, ...         % Reduce the number of partition boxes
    'TolXAbs', 1e-14, ...      % Relax absolute error tolerance
    'TolXRel', 1e-14, ...      % Relax relative error tolerance
    'NIT', 5, ...              % Reduce the number of global iterations
    'ND', 10 ...                % Reduce the number of local iterations
);

[X_inf, XS_inf] = verifynlssall(f, x0,options);

% Check if XS length exceeds 0
if length(XS_inf) > 0
    XS_inf
    error('XS length exceeds 0. Stopping execution.');
end

lower_bounds_X = inf(X_inf);
[~, idx] = sort(lower_bounds_X);
X_inf = X_inf(idx);

f = @(x)(func_left_hand_side(sup(s),x)); % Use singularity function if s < 0.05

[X_sup, XS_sup] = verifynlssall(f, x0,options);

% Check if XS length exceeds 0
if length(XS_sup) > 0
    XS_sup
    error('XS length exceeds 0. Stopping execution.');
end

lower_bounds_X = inf(X_sup);
[~, idx] = sort(lower_bounds_X);
X_sup = X_sup(idx);

zeros = I_hull(X_inf,X_sup); % Assign the computed intervals to zeros
zeros_for_each_s = arrayfun(@(i) [inf(zeros(i)), sup(zeros(i))], 1:length(zeros), 'UniformOutput', false);
intervals_flat = [zeros_for_each_s{:}]; % Flatten the intervals array

% Combine the index and the intervals into one row
result_for_each_s = [i, intervals_flat]

% Append the result to the CSV file
writematrix(result_for_each_s, csv_file, 'WriteMode', 'append');

%-----------------

for i=1:2^N1-1    
    s = hull(i-1,i)*bi_N1
    sigma = (1-s)^(1/intval(3));
    f = @(x)(func_left_hand_side_singularity(sup(sigma),x));
    % f = @(x)(func_left_hand_side(sup(s),x)); % Use singularity function if s < 0.05

    options = verifynlssallset( ...
        'Boxes', 2^10, ...         % Reduce the number of partition boxes
        'TolXAbs', 1e-14, ...      % Relax absolute error tolerance
        'TolXRel', 1e-14, ...      % Relax relative error tolerance
        'NIT', 5, ...              % Reduce the number of global iterations
        'ND', 10 ...                % Reduce the number of local iterations
    );

    [X, XS] = verifynlssall(f, x0,options);

    % Check if XS length exceeds 0
    if length(XS) > 0
        XS
        error('XS length exceeds 0. Stopping execution.');
    end

    lower_bounds_X = inf(X);
    [~, idx] = sort(lower_bounds_X);
    X = X(idx);

    X = error_left_hand_side(s,X);

    zeros = X; % Assign the computed intervals to zeros
    zeros_for_each_s = arrayfun(@(i) [inf(zeros(i)), sup(zeros(i))], 1:length(zeros), 'UniformOutput', false);
    intervals_flat = [zeros_for_each_s{:}]; % Flatten the intervals array

    % Combine the index and the intervals into one row
    result_for_each_s = [i, intervals_flat]

    % Append the result to the CSV file
    writematrix(result_for_each_s, csv_file, 'WriteMode', 'append');
end


for i=1:2^N2
    % s = hull(N1-2,N1-1)/N1+hull(1,N1)/N1*hull(i-1-N1,i-N1)/N2;
    s_inf = (2^N1-1)*bi_N1 + bi_N1*bi_N2*(i-1);
    s_sup = (2^N1-1)*bi_N1 + bi_N1*bi_N2*i;
    s = infsup(s_inf,s_sup);

    sigma = (1-s)^(1/intval(3));
    f = @(x)(func_left_hand_side_singularity(sup(sigma),x));

    options = verifynlssallset( ...
        'Boxes', 2^10, ...         % Reduce the number of partition boxes
        'TolXAbs', 1e-14, ...      % Relax absolute error tolerance
        'TolXRel', 1e-14, ...      % Relax relative error tolerance
        'NIT', 8, ...              % Reduce the number of global iterations
        'ND', 10 ...                % Reduce the number of local iterations
    );
    
    [X, XS] = verifynlssall(f, x0,options);

    % Check if XS length exceeds 0
    if length(XS) > 0
        XS
        error('XS length exceeds 0. Stopping execution.');
    end

    lower_bounds_X = inf(X);
    [~, idx] = sort(lower_bounds_X);
    X = X(idx);
    
    X = error_left_hand_side_singularity(sigma,X);

    zeros = X; % Assign the computed intervals to zeros
    zeros_for_each_s = arrayfun(@(i) [inf(zeros(i)), sup(zeros(i))], 1:length(zeros), 'UniformOutput', false);
    intervals_flat = [zeros_for_each_s{:}]; % Flatten the intervals array

    % Combine the index and the intervals into one row
    result_for_each_s = [i+2^N1-1, intervals_flat]

    % Append the result to the CSV file
    writematrix(result_for_each_s, csv_file, 'WriteMode', 'append');
end

