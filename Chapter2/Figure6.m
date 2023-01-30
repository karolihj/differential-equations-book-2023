% Compute the numerical solution of the FitzHugh-Nagumo model and make some
% parameter adjustments. Plot the upstroke and the action potential.
clear all

% Set up figure
figure('Units','centimeters', 'Position', [10 10 20 22], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [18, 20])

% Define time step
T = 5000;   % Simulation time
N = 5000;   % Number of time points
dt = T/N;   % Time step

parameter_names = {'a', 'c1', 'c2', 'b', 'd'};

for i=1:length(parameter_names)
    legends = {};

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

        % Update legends
        legends{j} = sprintf('%g\\times%s', adjust_values(j), ...
            parameter_names{i});

        % Set up arrays for saving the solutions
        t = (0:dt:T);  % Time vector
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
        
        % Time point of maximum upstroke velocity of the second half of the
        % solution
        [~, t_idx] = max((v(round(N/2)+1:end-round(N/4))-v(round(N/2):end-round(N/4)-1))/dt);
        t_idx = t_idx + round(N/2);
        
        % Plot the solution for a little while before the time of the
        % maximum upstroke velocity
        t_idx = t_idx - round(50/dt);
        t = t - t(t_idx);
        t = t(t_idx:end);
        v = v(t_idx:end);
        
        % Plot the upstroke
        subplot(5,3,(i-1)*3+1)
        plot(t, v, 'linewidth', 2)
        hold on
        set(gca, 'fontsize', 12)
        if i==1
            title('Upstroke')
        elseif i==5
            xlabel('t')
        end
        ylabel('v')
        xlim([20, 80])
        
        % Plot the action potential
        subplot(5,3,(i-1)*3+2)
        plot(t, v, 'linewidth', 2)
        hold on
        set(gca, 'fontsize', 12)
        if i==1
            title('Action potential')
        elseif i==5
            xlabel('t')
        end
        xlim([0, 530])
        
        % Plot legends
        subplot(5,3,(i-1)*3+3)
        plot(1, 1, 'linewidth', 2)
        hold on
        set(gca, 'fontsize', 12)
        axis off
        if j==length(adjust_values)
            legend(legends, 'location', 'westoutside', 'fontsize', 12)
        end
        xlim([0, 530])
    end
end

% Save figure
print('-dpdf', '../Figures/Ch2_Fig6.pdf')

