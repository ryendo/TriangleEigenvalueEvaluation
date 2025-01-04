format compact short infsup
close all



INTERVAL_MODE=1;


% Ensure INTLAB is initialized
% Add INTLAB to your MATLAB path and initialize it before using

% Define your function f
s = 0;  % Example value for s
s = I_infsup(0,0.0001)
f = @(x)(func_left_hand_side(s,x));
% f = @(x)(func_left_hand_side_singularity(s,x));

% Define the interval X excluding 0 to 0.15
X = infsup(0, 5);  % Adjust the lower bound to 0.15

% Number of subintervals
N = 100;

% Find intervals containing zeros
zero_intervals = find_zeros(f, X, N);

% Generate points for plotting from 0.15 to 10
x_plot = linspace(inf(X), sup(X), 1000);
y_plot = zeros(size(x_plot));

% Evaluate f at the plotting points
for i = 1:length(x_plot)
    xi = x_plot(i);
    yi = f(intval(xi));
    y_plot(i) = mid(yi);  % Use the midpoint for plotting
end

% Plot the function
figure;
plot(x_plot, y_plot, 'b-', 'LineWidth', 1.5);
hold on;

% Highlight intervals containing zeros in red
for i = 1:size(zero_intervals, 1)
    xi = zero_intervals(i);
    x_inf = inf(xi);
    x_sup = sup(xi);
    x_rect = [x_inf, x_sup, x_sup, x_inf];
    y_rect = ylim;
    y_rect = [y_rect(1), y_rect(1), y_rect(2), y_rect(2)];
    fill(x_rect, y_rect, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

% Add labels and title
xlabel('x');
ylabel('f(x)');
title('Function f(x) with Intervals Containing Zeros Highlighted in Red');
grid on;
xlim([0.15, sup(X)]);  % Set x-axis limits to exclude 0 to 0.15
hold off;