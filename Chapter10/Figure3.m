% Compute the solution of the bidomain model with membrane dynamics modeled
% by the parsimonious ventricular rabbit model and see how the cv depends
% on sigma_i and sigma_e
clear all

% Default values
sigma_i = 3;
sigma_e = 10;

% Varying values
sigma_i_values = [1, 2, 3, 4, 5];
sigma_e_values = [2, 5, 10, 15, 20];

% Perform simulations
cv_i = zeros(length(sigma_i_values), 1);
cv_e = zeros(length(sigma_e_values), 1);
for i=1:length(sigma_i_values)
    cv_i(i) = solve_system_and_compute_cv(sigma_i_values(i), sigma_e, ...
        2*(i-1)+1, 2*length(sigma_i_values));
    cv_e(i) = solve_system_and_compute_cv(sigma_i, sigma_e_values(i), ...
        2*(i-1)+2, 2*length(sigma_i_values));
end

% Set up figure and plot the conduction velocity
figure('Units','centimeters', 'Position', [1 1 22 7], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [19, 7]);

subplot(1,2,1)
plot(sigma_i_values, cv_i, '.-', 'linewidth', 2, 'Markersize', 20)
set(gca, 'fontsize', 12)
xlabel('\sigma_i (mS/cm)')
title('\sigma_i', 'fontsize', 18)
ylabel('CV (cm/s)')
ylim([20, 70])
xlim([min(sigma_i_values), max(sigma_i_values)])

subplot(1,2,2)
plot(sigma_e_values, cv_e, '.-', 'linewidth', 2, 'Markersize', 20)
set(gca, 'fontsize', 12)
xlabel('\sigma_e (mS/cm)')
title('\sigma_e', 'fontsize', 18)
ylabel('CV (cm/s)')
ylim([20, 70])
xlim([min(sigma_e_values), max(sigma_e_values)])

% Save the figure
print('-dpdf', '../Figures/Ch10_Fig3.pdf')