% Plot the analytical and the numerical solution of the problem
% y'(t) = y(t), y(0) = 1, for t between 0 and 1.
clear all

% Analytical solution
T = 1;
t_analytical = (0:0.01:T);
y_analytical = exp(t_analytical);

% Set up figure
figure('Units','centimeters', 'Position', [10 10 28 6], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [24, 6])

% Consider a number of different values of dt = T/N = 1/N
N_values = [5, 10, 100, 1000];

for i=1:length(N_values)
    
    % Numerical solution
    N = N_values(i);
    dt = T/N;
    t_numerical = (0:dt:T);
    y_numerical = (1+dt).^(0:N);
    
    % Plot the solutions
    subplot(1, length(N_values), i)
    plot(t_analytical, y_analytical, 'linewidth', 2)
    hold on
    plot(t_numerical, y_numerical, ':', 'linewidth', 2)
    set(gca, 'fontsize', 12)
    title(sprintf('\\Deltat = %g', dt))
    xlabel('t')
    if i==1
        ylabel('y')
    end
    if i==length(N_values)
        l = legend({'analytical', 'numerical'}, 'Location', 'northwest');
    end

end

print('-dpdf', '../Figures/Ch1_Fig1.pdf')