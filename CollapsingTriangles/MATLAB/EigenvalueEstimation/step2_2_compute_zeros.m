format long infsup
INTERVAL_MODE=1;

x0 = infsup(0,4);


%treat the first some intervals differently because dF/Dk vanishses at s=0

s = 0;
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

[X, XS] = verifynlssall(f, x0,options)
