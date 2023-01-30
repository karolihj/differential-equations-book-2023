% Compute the numerical solution of the FitzHugh-Nagumo model with an added
% diffusion term
clear all

% Define model constants
a = -0.12;
c1 = 0.175;
c2 = 0.03;
b = 0.011;
d = 0.55;
delta = 5e-5;

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


% Set up arrays for saving the solutions
v = zeros(M, N+1);
w = zeros(M, N+1);

% Define initial conditions
v(:,1) = 0;
w(:,1) = 0;

% Define non-zero v at the left part of the domain
v(1:5,1) = 0.26;

% Define diffusion matrix
D = spdiags((1-2*rho)*ones(M, 1), 0, M, M) + ...
    spdiags([0; 2*rho; rho*ones(M-2, 1)], 1, M, M) + ...
    spdiags([rho*ones(M-2, 1); 2*rho; 0], -1, M, M);

% Compute the numerical solution
for n=1:N
    v(:,n+1) = D*v(:,n) ...
        + dt*(c1*v(:,n).*(v(:,n)-a).*(1-v(:,n)) - c2*w(:,n));
    w(:,n+1) = w(:,n) + dt*(b*(v(:,n)-d*w(:,n)));
end

% Select time points to plot
t_plot = [0, 50, 100, 200, 300];
t_idx = round(t_plot/dt) + 1;

% Set up figure and plot the solution
figure('Units','centimeters', 'Position', [10 10 25 5], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [21, 5])

for i=1:length(t_plot)
    subplot(1, length(t_plot), i)
    plot(x, v(:, t_idx(i)), 'linewidth', 2)
    set(gca, 'fontsize', 12)
    title(sprintf('t = %g', t_plot(i)))
    xlabel('x')
    if i == 1
        ylabel('v')
        
    end
    ylim([-0.4, 1.4])
end

% Save figure
print('-dpdf', '../Figures/Ch6_Fig1.pdf')



