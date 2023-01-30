% Solve a non-linear reaction-diffusion equation using operator splitting 
% and compare to the explicit solution using very fine time stepping
clear all

% Print start of table
fprintf('dt      E         E/dt\n')

%% EXPLICIT SOLUTION USING VERY FINE TIME STEPPING
% Discretization parameters
L = 1;                % Spatial domain length
T = 1;                % Total simulation time
dx = 0.01;            % Distance between spatial grid points
x = (0:dx:1);         % x vector
dt = 1e-6;            % Time step
M = round(L/dx) + 1;  % Number of spatial grid points
N = round(T/dt);      % Number of time steps


% Set up diffusion matrix
D = (1/dx^2)*(spdiags(-2*ones(M-2, 1), 0, M-2, M-2) + ...
    spdiags(ones(M-2, 1), 1, M-2, M-2) + ...
    spdiags(ones(M-2, 1), -1, M-2, M-2));

% Set up identity matrix
I = speye(M-2);

% Solution vectors
u_fine = zeros(M, N+1);

% Initial conditions
u_fine(:, 1) = sin(pi*x); 

% Numerical time stepping (for all spatial points except the first and last)
for n=1:N
    u_fine(2:end-1, n+1) = (I + dt*D)*u_fine(2:end-1, n) - dt*u_fine(2:end-1, n).^2;
end


%% IMPLICIT OPERATOR SPLITTING SOLUTION
dt_values = [0.01, 0.005, 0.001, 0.0005, 0.0001];
for i=1:length(dt_values)
    dt = dt_values(i);
    N = round(T/dt);    % Number of time steps

    % Solution vector
    u = zeros(M, N+1);

    % Initial conditions
    u(:, 1) = sin(pi*x);

    % Numerical time stepping (for all spatial points except the first and last)
    for n=1:N
        % Step 1 (diffusion part)
        u(2:end-1, n+1) = (I - dt*D)\u(2:end-1, n);

        % Step 2 (reaction part)
        u(2:end-1, n+1) = (-1 + sqrt(1+4*dt*u(2:end-1, n+1)))/(2*dt);
    end

    % Compute maximum error
    error = max(abs(u(:, end)-u_fine(:, end)));
    
    % Print the error
    fprintf('%-7g %-9.3g %-9.2g\n', dt, error, error/dt);

end