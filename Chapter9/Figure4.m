% Solve the equations of the cable equation with membrane dynamics modeled
% by the Hodgkin-Huxley model using an operator splitting scheme and
% compute the conduction velocity for a few parameter changes
clear all

% Set up figure
figure('Units','centimeters', 'Position', [10 10 30 15], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [26, 15])

% Set up default parameters
w = 0.001;       % cm
sigma_i = 4;     % mS/cm
g_Na = 120;      % mS/cm^2
g_K = 36;        % mS/cm^2
g_L = 0.3;       % mS/cm^2

% Specify adjustments
titles = {'g_{Na}', 'g_{K}', 'g_L', 'w', '\sigma_i'};
units = {'mS/cm^2', 'mS/cm^2', 'mS/cm^2', 'cm', 'mS/cm'};

% Suitable adjustment factors for each parameter
factors = zeros(5, 11);
factors(1,:) = linspace(0.5, 1.2, 11);
factors(2,:) = linspace(0.8, 1.6, 11);
factors(3,:) = linspace(0.5, 1.5, 11);
factors(4,:) = linspace(0.5, 1.5, 11);
factors(5,:) = linspace(0.5, 1.5, 11);

for i=1:length(titles)

    % Run simulations
    CVs = zeros(length(factors), 1);
    for j=1:length(factors)
        p = [g_Na, g_K, g_L, w, sigma_i];
        p(i) = factors(i,j)*p(i);
        CVs(j) = solve_system_and_compute_cv(p(1), p(2), p(3), p(4), p(5), ...
            (i-1)*length(factors)+j, 5*length(factors));
    end

    % Plot the conduction velocity
    p0 = [g_Na, g_K, g_L, w, sigma_i];
    subplot(2,3,i)
    plot(factors(i,:)*p0(i), CVs, '.-', 'Linewidth', 2, 'Markersize', 20)
    set(gca, 'fontsize', 14)
    title(titles{i})
    xlabel(sprintf('%s (%s)', titles{i}, units{i}))
    if rem(i-1, 3) == 0
        ylabel('CV (cm/s)')
    end
    ylim([40, 90])
end

% Save figure
print('-dpdf', '../Figures/Ch9_Fig4.pdf')
