% Compute the numerical solution of the FitzHugh-Nagumo model with an added
% diffusion term - vary the value of c1
clear all

% Define model constants
a = -0.12;
c1 = 0.175;
c2 = 0.03;
b = 0.011;
d = 0.55;
delta = 5e-5;

c1_values = [0.1, 0.2, 0.4];

% Define discretization parameters
L = 1;                % Length of domain
dx = 0.01;            % Distance between spatial grid points
M = round(L/dx) + 1;  % Number of spatial grid points
x = (0:dx:L);         % x vector
dt = 0.005;           % Time step
T = 300;              % Simulation time
N = round(T/dt);      % Number of time points 
t = (0:dt:T);         % Time vector
rho = delta*dt/dx^2;

% Define diffusion matrix
D = spdiags((1-2*rho)*ones(M, 1), 0, M, M) + ...
    spdiags([0; 2*rho; rho*ones(M-2, 1)], 1, M, M) + ...
    spdiags([rho*ones(M-2, 1); 2*rho; 0], -1, M, M);

% Set up figure
figure('Units','centimeters', 'Position', [10 10 18 24], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [17, 23])

% Select time points to plot
t_plot = [0, 50, 100, 150, 200];
t_idx = round(t_plot/dt) + 1;

for k=1:length(c1_values)
    c1 = c1_values(k);

    % Set up arrays for saving the solutions
    v = zeros(M, N+1);
    w = zeros(M, N+1);

    % Define initial conditions
    v(:,1) = 0;
    w(:,1) = 0;
    
    % Define non-zero v at the left part of the domain
    v(1:5,1) = 0.26;

    % Compute the numerical solution
    for n=1:N
        v(:,n+1) = D*v(:,n) ...
            + dt*(c1*v(:,n).*(v(:,n)-a).*(1-v(:,n)) - c2*w(:,n));
        w(:,n+1) = w(:,n) + dt*(b*(v(:,n)-d*w(:,n)));
    end

    % Plot the numerical solution
    for i=1:length(t_plot)
        subplot(length(t_plot), length(c1_values), ...
            (i-1)*length(c1_values) + k)
        plot(x, v(:, t_idx(i)), 'linewidth', 2)
        set(gca, 'fontsize', 14)
        text(0.98, 0.88, sprintf('t = %g', t_plot(i)), 'FontSize', 14, ...
            'Units','normalized', 'HorizontalAlignment','right')
        if i==1
            title(sprintf('c_1 = %g', c1_values(k)))
        end
        if k==1
            ylabel('v')
        end
        if i == length(t_plot)
            xlabel('x')
        end
        ylim([-0.4, 1.4])
    end
end

% Save figure
print('-dpdf', '../Figures/Ch6_Fig3.pdf')