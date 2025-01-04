%% [Numerical version] 基底: (x - aL) * (aR - x) * a * (x - c_i)^2
% clear; clc;  % 必要に応じて

%====================================
% 1) パラメータ設定
%====================================
N = 16;                 % 例: 基底関数の個数
a_val = 1;           % パラメータ a (double)
t_val = 2*tan(pi/18);  % tの値 (double)
% t_val = 0.1;  % tの値 (double)
s_val = 0.5;           % sの値 (double)

%--------------------------------------------
% 端点 a_left, a_right, および分割点 x_split
%--------------------------------------------
a_left_num  = - t_val^(-2/3)*(1 + s_val);
a_right_num =   t_val^(-2/3)*(1 - s_val);
x_split_num = - s_val / (t_val^(2/3));

fprintf('N=%d\n', N);
fprintf('Interval: [%.6f, %.6f], split = %.6f\n', ...
    a_left_num, a_right_num, x_split_num);

%====================================
% 2) 基底関数 phi_i(x), dphi_i(x) の定義
%====================================
% c_i = a_left + (i-1)*((a_right - a_left)/(N-1))
c_i = @(i) a_left_num + (i-1)*((a_right_num - a_left_num)/(N-1));

% Define the basis function phi_i(x) using sine functions
phi_i = @(x, i) sin(i * pi * (x - a_left_num) / (a_right_num - a_left_num));

% Define the derivative of the basis function dphi_i/dx
dphi_i = @(x, i) (i * pi / (a_right_num - a_left_num)) * cos(i * pi * (x - a_left_num) / (a_right_num - a_left_num));


% phi_i = @(x, i)(x - a_left_num).*(a_right_num - x).* a_val .* (x - c_i(i)).^2;
% phi_i = @(x, i) ...
%     (x - a_left_num) .* (a_right_num - x) .* a_val .* (x - c_i(i)).^5;
% phi_i = @(x, i) ...
%     (x - a_left_num) .* (a_right_num - x) .* exp(-a_val .* (x - c_i(i)).^2);
% phi_i = @(x, i) (x >= a_left_num + (i-1)*(a_right_num - a_left_num)/N) & ...
%                 (x < a_left_num + i*(a_right_num - a_left_num)/N);

% taylor_order = 5;
% syms z
% phi_i = @(x, i) ...
%     (x - a_left_num) .* (a_right_num - x) .* ...
%     polyval(fliplr(double(coeffs(taylor(exp(-a_val .* (z - c_i(i)).^2), z, 'Order', taylor_order)))), x - c_i(i));

% % 区分的線形基底関数
% phi_i = @(x, i) max(0, min((x - c_i(i-1))/(c_i(i) - c_i(i-1)), ...
%                            (c_i(i+1) - x)/(c_i(i+1) - c_i(i))));
% 
% % 微分 dphi_i/dx は解析的に実装
% dphi_i = @(x, i) ...
%     (x > c_i(i-1) & x <= c_i(i)) .* (1/(c_i(i) - c_i(i-1))) + ...
%     (x > c_i(i) & x <= c_i(i+1)) .* (-1/(c_i(i+1) - c_i(i)));

% Define the basis function phi_i(x) as a polynomial
phi_i = @(x, i) (x - a_left_num) .* (a_right_num - x) .* a_val .* x.^(i-1);

% Define the derivative of the basis function dphi_i/dx
dphi_i = @(x, i) a_val * ( ...
    (a_right_num - x) .* (i-1) .* x.^(i-2) - (x - a_left_num) .* (i-1) .* x.^(i-2) ...
);

% 微分 dphi_i/dx は数値微分(中央差分)で定義 (実用上は解析的微分のほうが望ましい)
% h = 1e-12;
% dphi_i = @(x, i) (phi_i(x+h, i) - phi_i(x-h, i))./(2*h);

%====================================
% 3) V_minus, V_plus の定義
%====================================
% 定義に登場する h_p(x), h_m(x) など
% Define V_minus, V_plus
Y_fun = @(x) t_val^(2/3).*x + s_val;

h_p  = @(x) ( -t_val / (1 - s_val) ) .* ( Y_fun(x) - 1 );
dh_p = @(x) ( -t_val / (1 - s_val) );

h_m  = @(x) ( t_val / (1 + s_val) ) .* ( Y_fun(x) + 1 );
dh_m = @(x) ( t_val / (1 + s_val) );

V_plus = @(x) t_val^(4/3) .* ( ...
    ( pi^2 ./ (h_p(x).^2) ) ...
  + ( (3 + 4*pi^2) .* (dh_p(x).^2 ) ) ./ ( 12 .* (h_p(x).^2) ) ...
  - ( pi^2 / t_val^2 ) ...
);

V_minus = @(x) t_val^(4/3) .* ( ...
    ( pi^2 ./ (h_m(x).^2) ) ...
  + ( (3 + 4*pi^2) .* (dh_m(x).^2 ) ) ./ ( 12 .* (h_m(x).^2) ) ...
  - ( pi^2 / t_val^2 ) ...
);

%====================================
% 4) 行列 K, M の組み立て (すべて数値積分)
%====================================
K_mat = zeros(N,N);
M_mat = zeros(N,N);

% 端点
ep = 0.01;
xL = a_left_num+ep;
xR = a_right_num-ep;
xC = x_split_num;

% 数値積分を行う補助関数
numIntSplit = @(f) integral(f, xL, xC) + integral(f, xC, xR);

for iInd = 1:N
    for jInd = iInd:N  % 対称なので jInd=iIndから

        % K_ij:
        %   (dphi_i * dphi_j) + V_- * (phi_i * phi_j) on [xL,xC]
        %   (dphi_i * dphi_j) + V_+ * (phi_i * phi_j) on [xC,xR]
        integrandK_left  = @(x) dphi_i(x,iInd).*dphi_i(x,jInd) ...
                                + V_minus(x).*phi_i(x,iInd).*phi_i(x,jInd);
        integrandK_right = @(x) dphi_i(x,iInd).*dphi_i(x,jInd) ...
                                + V_plus(x).*phi_i(x,iInd).*phi_i(x,jInd);

        K_val = integral(integrandK_left, xL, xC) ...
              + integral(integrandK_right, xC, xR);

        % M_ij:
        %   phi_i * phi_j を区間わけなしで2回足すだけ(中身同じなので V_±関係なし)
        integrandM = @(x) phi_i(x,iInd).*phi_i(x,jInd);
        M_val = integral(integrandM, xL, xC) ...
              + integral(integrandM, xC, xR);

        K_mat(iInd,jInd) = K_val;
        K_mat(jInd,iInd) = K_val;  % 対称行列

        M_mat(iInd,jInd) = M_val;
        M_mat(jInd,iInd) = M_val;
    end
end

%====================================
% 5) 固有値問題を解く
%====================================
fprintf('Solving generalized eigenvalue problem K * v = lambda * M * v ...\n');
[V_approx, D_approx] = eig(K_mat, M_mat);
lambda_approx = diag(D_approx);

[lambda_sorted, idx] = sort(lambda_approx);
V_sorted = V_approx(:, idx);

% 下位2つ(あるいは小さい方から2つ)を表示
disp('Approx eigenvalues:');
disp(lambda_sorted(1:2));

%====================================
% 6) 固有関数(近似)を作図
%====================================
k_plot = 1;  % 2番目の固有値に対応する固有ベクトルを例に

x_plot = linspace(xL, xR, 300);
u_plot = zeros(size(x_plot));

for iInd = 1:N
    coeff = V_sorted(iInd, k_plot);
    u_plot = u_plot + coeff * phi_i(x_plot, iInd);
end

% L2ノルムで正規化（任意）
norm_u = sqrt(trapz(x_plot, u_plot.^2));
if norm_u > 1e-14
    u_plot = u_plot / norm_u;
end

figure;
plot(x_plot, u_plot, 'b-', 'LineWidth',2);
xlabel('x'); 
ylabel(sprintf('u_{%d}(x)', k_plot));
title(sprintf('Approx. Eigenfunction %d (N=%d)', k_plot, N));
grid on;
