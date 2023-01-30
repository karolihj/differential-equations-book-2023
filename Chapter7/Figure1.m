% Solve a non-linear reaction-diffusion equation using operator splitting
clear all

% Discretization parameters
L = 1;                % Length of domain
T = 1;                % Total simulation time  
dx = 0.01;            % Distance between spatial grid points
x = (0:dx:1);         % xvector
dt = 0.001;           % Time step
M = round(L/dx) + 1;  % Number of spatial grid points
N = round(T/dt);      % Number of time steps


% Set up diffusion matrix
D = (1/dx^2)*(spdiags(-2*ones(M-2, 1), 0, M-2, M-2) + ...
    spdiags(ones(M-2, 1), 1, M-2, M-2) + ...
    spdiags(ones(M-2, 1), -1, M-2, M-2));

% Set up identity matrix
I = speye(M-2);

% Solution vector
u = zeros(M, N+1);

% Initial conditions
u(:, 1) = sin(pi*x); 

% Numerical time stepping (for all spatial points except the first
% and last, which are known from the boundary conditions)
for n=1:N

    % Step 1 (diffusion part)
    u(2:end-1, n+1) = (I - dt*D)\u(2:end-1, n);

    % Step 2 (reaction part)
    u(2:end-1, n+1) = (-1 + sqrt(1+4*dt*u(2:end-1, n+1)))/(2*dt);

end


% Set up figure and plot the solution at time T
figure('Units','centimeters', 'Position', [10 10 15 10], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [14, 10])

plot(x, u(:, end), 'linewidth', 3)
set(gca, 'fontsize', 14)
ylabel('u')
xlabel('x')

% Save figure
print('-dpdf', '../Figures/Ch7_Fig1.pdf')