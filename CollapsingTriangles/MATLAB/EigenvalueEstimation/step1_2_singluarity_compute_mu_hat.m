format compact short infsup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Symbolic definitions with N as a numeric integer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We do NOT do "clear; clc;" per your requirement, so watch for conflicts.
clear

% ---------------- NEW: set N as a numeric integer ----------------
N = 15;  % for example

% Symbolic variables: i, j, plus t, s (but NOT N)
syms i j t sigma a real
assume(i, 'integer');
assume(j, 'integer');
assume(t > 0);
assume(abs(sigma) < 1);

% We will store them in one file, e.g. "K_M_sym_withN.mat"
symFilename = 'upper_bound_matrix/K_M_sym_withN_5_xi_NN_mfacnn.mat';

if exist(symFilename, 'file')
    fprintf('Loading symbolic integrals from file: %s\n', symFilename);
    load(symFilename);
% if false
%     fprintf('Loading symbolic integrals from file: %s\n', symFilename);
%     load(symFilename);
else
    fprintf('No existing symbolic integrals file. Building them now...\n');

    % --- Define a_left, a_right, x_split in symbolic form ---
    a_left_sym  = -t^(-sym(2)/sym(3))*(1 + s);
    a_right_sym =  t^(-sym(2)/sym(3))*(1 - s);
    x_split_sym = -s/(t^(sym(2)/sym(3)));   % solves t^(2/3)*x + s = 0

    % Symbolic x for integration
    syms x real

    % Define V_minus, V_plus
    % Y_sym = t^(2/sym(3))*x + s;
    % h_p_sym       = @(Y) -t/(1 - s)*(Y - 1);
    % h_m_sym       = @(Y)  t/(1 + s)*(Y + 1);
    % h_prime_p_sym = @(Y) -t/(1 - s);
    % h_prime_m_sym = @(Y)  t/(1 + s);
    % 
    % V_plus_sym = ...
    %      t^(4/sym(3))*((sym(pi)^2 / (h_p_sym(Y_sym)^2)) ...
    %    + ((3 + 4*pi^2)*(h_prime_p_sym(Y_sym)^2)) / (12*(h_p_sym(Y_sym)^2)) ...
    %    - (sym(pi)^2 / t^2));
    % 
    % V_minus_sym = ...
    %      t^(4/sym(3))*((sym(pi)^2 / (h_m_sym(Y_sym)^2)) ...
    %    + ((3 + 4*pi^2)*(h_prime_m_sym(Y_sym)^2)) / (12*(h_m_sym(Y_sym)^2)) ...
    %    - (sym(pi)^2 / t^2));
    
    
    pi_sym = sym(pi);
    s = 1-sigma^3;
    % Define the numerator and denominator for the case x > 0
    numerator_pos = -24*pi_sym^2*(s-1)*x - 12*pi_sym^2*t^(2/3)*x^2 + (3 + 4*pi_sym^2)*t^(4/3);
    denominator_pos = 12*((s-1)^2 + 2*(s-1)*t^(2/3)*x + t^(4/3)*x^2);
    V_plus_sym = numerator_pos/denominator_pos;
    
    % Define the numerator and denominator for the case x <= 0
    numerator_neg = -24*pi_sym^2*s*x - 12*pi_sym^2*t^(2/3)*x^2 - 24*pi_sym^2*x + (3 + 4*pi_sym^2)*t^(4/3);
    denominator_neg = 12*((s+1)^2 + 2*(s+1)*t^(2/3)*x + t^(4/3)*x^2);
    V_minus_sym = numerator_neg/denominator_neg;

    % --- Define c_sym: the "node" for the i-th basis function
    c_sym = symfun( ...
        a_left_sym + (i - 1)*((a_right_sym - a_left_sym)/(N - 1)), ...
        i ...
    );

    % ==================================================================
    %     *****  Example polynomial basis with (x-aL)*(aR-x)*(x-c_i)^2  *****
    % ==================================================================
    % phi_sym = symfun( ...
    % nchoosek(N-1, i-1) * ((x - a_left_sym)/(a_right_sym - a_left_sym))^(i-1) * ...
    % (1 - (x - a_left_sym)/(a_right_sym - a_left_sym))^(N-i), ...
    % [x, i] ...
    % );
    % phi_sym = symfun( ...
    %     (x - a_left_sym) * (a_right_sym - x) * a*(x - c_sym(i)), ...
    %     [x, i] ...
    % );
    phi_sym = symfun( ...
        (x - a_left_sym) * (a_right_sym - x) * exp(-(x - c_sym(idx1))^2), ...
        [x, idx1] ...  
    );
    % phi_sym = symfun( ...
    % piecewise( ...
    %     x >= c_sym(i-1) & x < c_sym(i), ...
    %     (x - c_sym(i-1)) / (c_sym(i) - c_sym(i-1)), ...
    %     x >= c_sym(i) & x < c_sym(i+1), ...
    %     (c_sym(i+1) - x) / (c_sym(i+1) - c_sym(i)), ...
    %     0 ...
    % ), ...
    % [x, i] ...
    % );
    % phi_sym = symfun( (x - a_left_sym) .* (a_right_sym - x) .* x.^(i-1),[x,i]);
    % phi_sym = symfun((x - a_left_sym).^(i-1) .* (a_right_sym - x).^(N-i+1),[x,i]);

    % phi_sym = symfun( ...
    % (x - a_left_sym) .* (a_right_sym - x) .* ...
    % legendreP(i, 2*(x - a_left_sym)/(a_right_sym - a_left_sym) - 1), ...
    % [x, i] ...
    % );
    
    % T_sym = @(n, x) chebyshevT(n, x);
    % scaling_transform = @(x, a_left, a_right) (2*x - (a_left + a_right)) / (a_right - a_left);
    % phi_sym = symfun( ...
    % T_sym(i-1, scaling_transform(x, a_left_sym, a_right_sym)), ...
    % [x, i] ...
    % );

    % %%%%%%%%%%%%
    % syms z
    % taylor_order = 12; % Set the Taylor expansion order
    % gaussian_taylor_sym = taylor(exp(-a*(z - c_sym(i))^2), z, 'Order', taylor_order);
    % 
    % % Define phi_sym using the Taylor-expanded Gaussian
    % phi_sym = symfun( ...
    %     (x - a_left_sym) * (a_right_sym - x) * subs(gaussian_taylor_sym, z, x), ...
    %     [x, i] ...
    % );
    % %%%%%%%%%%%%
    
    dphi_sym = symfun( diff(phi_sym(x,idx1), x), [x, idx1]);

    % For i, j in {1..N}, define integrals for K, M
    K_mat_sym = sym(zeros(N,N));
    M_mat_sym = sym(zeros(N,N));

    phi_i_expr  = simplify(phi_sym(x, iInd));
    phi_j_expr  = simplify(phi_sym(x, jInd));
    dphi_i_expr = simplify(dphi_sym(x, iInd));
    dphi_j_expr = simplify(dphi_sym(x, jInd));

    % phi_i_expr  = phi_sym(x, iInd);
    % phi_j_expr  = phi_sym(x, jInd);
    % dphi_i_expr = dphi_sym(x, iInd);
    % dphi_j_expr = dphi_sym(x, jInd);

    integrand_K_left  = simplify(dphi_i_expr .* dphi_j_expr + V_minus_sym.*phi_i_expr.*phi_j_expr);
    integrand_K_right = simplify(dphi_i_expr .* dphi_j_expr + V_plus_sym .*phi_i_expr.*phi_j_expr);
    integrand_M_left  = simplify(phi_i_expr.*phi_j_expr);
    integrand_M_right = simplify(phi_i_expr.*phi_j_expr);

    % integrand_K_left  = dphi_i_expr .* dphi_j_expr + V_minus_sym.*phi_i_expr.*phi_j_expr;
    % integrand_K_right = dphi_i_expr .* dphi_j_expr + V_plus_sym .*phi_i_expr.*phi_j_expr;
    % integrand_M_left  = phi_i_expr.*phi_j_expr;
    % integrand_M_right = phi_i_expr.*phi_j_expr;

    K_left  = int(integrand_K_left,  x, a_left_sym,  x_split_sym);
    K_right = int(integrand_K_right, x, x_split_sym, a_right_sym);
    M_left  = int(integrand_M_left,  x, a_left_sym,  x_split_sym);
    M_right = int(integrand_M_right, x, x_split_sym, a_right_sym);

    valK = simplify(K_left + K_right);
    valM = simplify(M_left + M_right);

    for iInd_val = 1:N
        for jInd_val = iInd_val:N  % exploit symmetry

            % valK = rewrite(valK, 'erf');
            % valM = rewrite(valM, 'erf');

            % valK = K_left + K_right;
            % valM = M_left + M_right;

            K_mat_sym(iInd_val,jInd_val) = subs(valK, [iInd, jInd], [iInd_val, jInd_val]);
            K_mat_sym(jInd_val,iInd_val) = subs(valK, [iInd, jInd], [iInd_val, jInd_val]);  % symmetric
            M_mat_sym(iInd_val,jInd_val) = subs(valM, [iInd, jInd], [iInd_val, jInd_val]);
            M_mat_sym(jInd_val,iInd_val) = subs(valM, [iInd, jInd], [iInd_val, jInd_val]);
        end
    end

    % [commonFactor, K_mat_sym, M_mat_sym] = normalizeMatrices(K_mat_sym,M_mat_sym);
    % [commonFactor, K_mat_sym, M_mat_sym] = normalizeMatrices(K_mat_sym,M_mat_sym);

    % factor = K_mat_sym(N,N);
    % K_mat_sym=K_mat_sym*t^23;
    % M_mat_sym=M_mat_sym*t^23;

    K_sym = matlabFunction(K_mat_sym,'Vars',[a,s,t]);
    M_sym = matlabFunction(M_mat_sym,'Vars',[a,s,t]);

    fprintf('Saving symbolic integrals to file: %s\n', symFilename);
    save(symFilename);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) Build numeric interval matrices, solve eigenvalue problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nWe have N=%d from the code.\n', N);

a_val = intval(1);
t_val  = intval(2*tan(pi/18));
s_val  = intval(0.98);
sigma_val = (1-s_val)^(1/intval(3));

fprintf('\nBuilding K_mat, M_mat for N=%d\n', N);

% % Evaluate numeric endpoints for info
a_left_num_handle   = matlabFunction(-t^(-sym(2)/sym(3))*(1 + s));
a_right_num_handle  = matlabFunction(t^(-sym(2)/sym(3))*(1 - s));
x_split_num_handle  = matlabFunction(-s/(t^(sym(2)/sym(3))));
% 
a_left_num  = a_left_num_handle(s_val,t_val);
a_right_num = a_right_num_handle(s_val,t_val);
x_split_num = x_split_num_handle(s_val,t_val);
% 
% fprintf('Interval: [%.6f, %.6f], split=%.6f\n', mid(a_left_num), mid(a_right_num), mid(x_split_num));

% Build midpoint matrices for approximate solve
K_mat = I_zeros(N,N);
M_mat = I_zeros(N,N);

K_val_sym = K_sym(a_val ,sigma_val, t_val);  % pass midpoints as double
M_val_sym = M_sym(a_val ,sigma_val, t_val);


for ii = 1:N
    for jj = 1:N
        K_mat(ii,jj) = K_val_sym(ii,jj);
        M_mat(ii,jj) = M_val_sym(ii,jj);
    end
end

% (If you want interval matrices, you can do a small bounding around K_mat_mid, M_mat_mid.)

% Solve standard numeric problem at the midpoint:
veigs(K_mat, M_mat,'sm',2)
[V_approx, D_approx] = eig(mid(K_mat), mid(M_mat));

lambda_approx = diag(D_approx);
[lambda_sorted, idx] = sort(lambda_approx);
V_sorted = V_approx(:, idx);

disp('Approx eigenvalues (midpoint):');
disp(lambda_sorted(1:3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3) Plot an approximate eigenfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Let's pick an eigenfunction to plot
k_plot = 2;  % for example, the 2nd eigenvalue

% Make a 1D grid
xL = double(mid(a_left_num));
xR = double(mid(a_right_num));
xa = double(I_mid(a_val));
x_plot = linspace(xL, xR, 300);

% We'll reconstruct
u_plot = zeros(size(x_plot));

% We *must* define a numeric version of our basis function:
phi_i_fun = matlabFunction(phi_sym(x,i));
% where c_i_num = a_left_num + (i-1)*((a_right_num - a_left_num)/(N-1)).

for iInd = 1:N
    coeff = V_sorted(iInd, k_plot);
    u_plot = u_plot + coeff * phi_i_fun(iInd,sigma_val,t_val, x_plot);
end

% Normalize if desired
L2_norm = sqrt(trapz(x_plot, u_plot.^2));
if L2_norm > 1e-14
    u_plot = u_plot / L2_norm;
end

figure;
plot(x_plot, u_plot, 'b-', 'LineWidth',2);
xlabel('x');
ylabel(sprintf('u_{%d}(x)', k_plot));
title(sprintf('Approx. Eigenfunction %d (N=%d)', k_plot, N));
grid on;