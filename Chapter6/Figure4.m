% Compute the numerical solution of the FitzHugh-Nagumo model with an added
% diffusion term and compute the conduction velocity 
% - vary the value of c1 and delta
clear all

% Default parameter values
c1 = 0.175;
delta = 5e-5;

% Varying parameter values
c1_values = 0.1:0.05:0.4;
delta_values = 1e-5:2e-5:20e-5;

% Perform simulations and compute the conduction velocity
CV_c1 = zeros(length(c1_values), 1);
CV_delta = zeros(length(delta_values), 1);

for i=1:length(c1_values)
    CV_c1(i) = solve_system_and_compute_cv(c1_values(i), delta);
end

for i=1:length(delta_values)
    CV_delta(i) = solve_system_and_compute_cv(c1, delta_values(i));
end


% Set up figure and plot the conduction velocity
figure('Units','centimeters', 'Position', [10 10 26 9], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [25, 9])

subplot(1,2,1)
plot(c1_values, CV_c1, '.-', 'linewidth', 2, 'Markersize', 24)
set(gca, 'fontsize', 14)
xlabel('c_1')
title('c_1')
ylabel('CV')

subplot(1,2,2)
plot(delta_values, CV_delta, '.-',  'linewidth', 2, 'Markersize', 24)
set(gca, 'fontsize', 14)
title('\delta')
xlabel('\delta')

% Save figure
print('-dpdf', '../Figures/Ch6_Fig4.pdf')