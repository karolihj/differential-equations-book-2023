% Solve the diffusion equation using an explicit numerical scheme 
% with a large dt
clear all

% Discretization parameters
T = 0.5;       % Total simulation time
M = 51;        % Number of spatial points
dx = 1/(M-1);  % Distance between spatial grid points
N = 500;       % Number of time points
dt = T/N;      % Time step
rho = dt/dx^2;

% Solution vector
u = zeros(M, N+1);

% Initial conditions
u(1:round(M/2), 1) = 1;
u(round(M/2)+1:end, 1) = 0;

% Numerical scheme
for n=1:N
    for j=1:M
        if j==1 % Left boundary
            u(j, n+1) = (1-2*rho)*u(j,n) + 2*rho*u(j+1,n);
        elseif j==M % Right boundary
            u(j, n+1) = 2*rho*u(j-1,n) + (1-2*rho)*u(j,n);
        else % Inner point
            u(j, n+1) = rho*u(j-1,n) + (1-2*rho)*u(j,n) + rho*u(j+1,n);
        end
    end
end

% Select time points to plot
x = (0:dx:1);
t_plot = [0, 0.001, 0.01, 0.05, 0.1];
t_idx = round(t_plot/dt) + 1;

% Set up figure and plot the solution
figure('Units','centimeters', 'Position', [10 10 10 18], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [10, 17])

for i=1:length(t_plot)
    subplot(length(t_plot), 1, i)
    plot(x, u(:, t_idx(i)), 'linewidth', 2)
    set(gca, 'fontsize', 12)
    title(sprintf('t = %g', t_plot(i)))
    ylabel('u')
    if i == length(t_plot)
        xlabel('x')
    end
end

% Save figure
print('-dpdf', '../Figures/Ch3_Fig3.pdf')