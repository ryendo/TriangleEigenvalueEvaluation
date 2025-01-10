
% Define symbolic variables
syms t;
K = matlabFunction(t^(sym(4)/sym(3)));
% Substitute t with intval('0.1')
t_val = intval('0.1');
K_eval = K(t_val);


% K_original = @(t)t.^(4.0./3.0); % 元の関数ハンドル
% K_modified = @(t)K_original(intval(t)); % t と内部定数を自動的に intval に変換
% t_val = intval(0.1);
% K_modified(t_val)