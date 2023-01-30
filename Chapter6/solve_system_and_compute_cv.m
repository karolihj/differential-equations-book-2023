function CV = solve_system_and_compute_cv(c1, delta)
%CV = solve_system_and_compute_cv(c1, delta)

% Define model constants
a = -0.12;
c2 = 0.03;
b = 0.011;
d = 0.55;

% Define discretization parameters
L = 1;                % Length of domain
dx = 0.01;            % Distance between spatial grid points
M = round(L/dx) + 1;  % Number of spatial grid points
x = (0:dx:L);         % x vector
dt = 0.005;           % Time step
T = 1000;              % Simulation time
N = round(T/dt);      % Number of time points 
t = (0:dt:T);         % Time vector
rho = delta*dt/dx^2;

% Define diffusion matrix
D = spdiags((1-2*rho)*ones(M, 1), 0, M, M) + ...
    spdiags([0; 2*rho; rho*ones(M-2, 1)], 1, M, M) + ...
    spdiags([rho*ones(M-2, 1); 2*rho; 0], -1, M, M);

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
        v(:,n+1) = D*v(:, n) ...
            + dt*(c1*v(:,n).*(v(:,n)-a).*(1-v(:,n)) - c2*w(:,n));
        w(:,n+1) = w(:,n) + dt*(b*(v(:,n)-d*w(:,n)));
end

% Compute CV
x_start = 0.5;
x_end = 0.7;
v_th = 0.5;
CV = compute_cv(v, v_th, t, x, x_start, x_end);

end