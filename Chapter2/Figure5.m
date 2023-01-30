% Compute the numerical solution of the FitzHugh-Nagumo model and make some
% parameter adjustments
clear all

% Set up figure
figure('Units','centimeters', 'Position', [10 10 30 22], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [26, 20])

% Define time step
T = 5000;     % Simulation time
N = 5000;     % Number of time points
dt = T/N;     % Time step
t = (0:dt:T); % Time vector

parameter_names = {'a', 'c1', 'c2', 'b', 'd'};

for i=1:length(parameter_names)
   
    % Select adjustment factors suitable for each parameter
    if i == 2 || i == 5
        adjust_values = 0.6:0.2:1;
    elseif i == 3
        adjust_values = 1:0.2:1.4;
    else
        adjust_values = 0.6:0.4:1.4;
    end
    
    for j=1:length(adjust_values)
        
        % Define model constants
        a = -0.12;
        c1 = 0.175;
        c2 = 0.03;
        b = 0.011;
        d = 0.55;
        
        % Adjust parameters
        p = [a, c1, c2, b, d];
        p(i) = p(i)*adjust_values(j);
        a = p(1); c1 = p(2); c2 = p(3);  b = p(4); d = p(5);

        % Set up arrays for saving the solutions
        v = zeros(N+1, 1);
        w = zeros(N+1, 1);

        % Define initial conditions
        v(1) = 0.26;
        w(1) = 0;

        % Compute the numerical solution
        for n=1:N
            v(n+1) = v(n) + dt*(c1*v(n)*(v(n)-a)*(1-v(n)) - c2*w(n));
            w(n+1) = w(n) + dt*(b*(v(n)-d*w(n)));
        end
        
        % Plot the solution
        subplot(5,3,(i-1)*3+j)
        plot(t, v, 'linewidth', 2)
        set(gca, 'fontsize', 12)
        if i==5
            xlabel('t')
        end
        if j==1
            ylabel('v')
        end
        xlim([0, 5000])
        title(sprintf('%g\\times%s', adjust_values(j), ...
            parameter_names{i}))
    end
end

% Save figure
print('-dpdf', '../Figures/Ch2_Fig5.pdf')

