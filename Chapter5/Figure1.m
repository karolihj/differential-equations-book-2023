% Solve a linear reaction-diffusion equation using an explicit, an implicit 
% and a midpoint scheme
clear all

% Define lambda
lambda = 2*pi^2;

% Discretization parameters
T = 1;         % Total simulation time
M = 101;       % Number of spatial grid points
dx = 1/(M-1);  % Distance between spatial grid points
x = (0:dx:1);  % x-values
N = 1000;      % Number of time steps
dt = T/N;      % Time step

% Set up diffusion matrix
D = (1/dx^2)*(spdiags(-2*ones(M-2, 1), 0, M-2, M-2) + ...
    spdiags(ones(M-2, 1), 1, M-2, M-2) + ...
    spdiags(ones(M-2, 1), -1, M-2, M-2));

% Set up identity matrix
I = speye(M-2);

% Solution vectors
u_e = zeros(M, N+1);
u_i = zeros(M, N+1);
u_m = zeros(M, N+1);

% Initial conditions
u_e(:, 1) = sin(pi*x); 
u_i(:, 1) = sin(pi*x); 
u_m(:, 1) = sin(pi*x); 

% Numerical time stepping (for all spatial points except for the first 
% and last, which are known from the boundary conditions)
for n=1:N

    % Explicit scheme
    u_e(2:end-1, n+1) = ((1 + lambda*dt)*I + dt*D)*u_e(2:end-1, n);

    % Implicit scheme
    u_i(2:end-1, n+1) = ((1 - lambda*dt)*I - dt*D)\u_i(2:end-1, n);

    % Midpoint scheme
    u_m(2:end-1, n+1) = ((1 - lambda*dt/2)*I - (dt/2)*D)\(((1 + lambda*dt/2)*I + (dt/2)*D)*u_m(2:end-1, n));

end

% Analytical solution
u_a = exp(pi^2*T)*sin(pi*x);


% Set up figure and plot the solution at time T
figure('Units','centimeters', 'Position', [10 10 20 15], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [19, 15])
subplot(2,1,1)
plot(x, u_a, 'linewidth', 2)
hold on
plot(x, u_i(:, end), '--', 'linewidth', 2)
plot(x, u_m(:, end), ':', 'linewidth', 2)
set(gca, 'fontsize', 14)
ylabel('u')
xlabel('x')
legend({'Analytical solution', 'Implicit scheme', 'Midpoint scheme'}, ...
    'Location', 'south')

% Plot the explicit solution
t_idx = round(0.1/dt)+1;
subplot(2,1,2)
plot(x, u_e(:, t_idx), 'linewidth', 2, 'color', [0.4940 0.1840 0.5560])
set(gca, 'fontsize', 14)
ylabel('u')
xlabel('x')
legend({'Explicit scheme'}, 'Location', 'northwest')

% Save figure
print('-dpdf', '../Figures/Ch5_Fig1.pdf')