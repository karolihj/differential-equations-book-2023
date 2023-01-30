% Compute the error of the temporal part of the numerical solution of the
% diffusion equation at T = 1.
clear all

% Analytical solution
T = 1;
u_analytical = exp(-pi^2*T);

% Print start of table
fprintf('N        dt      E         E/dt\n')

% Consider a number of different values of dt = T/N = 1/N
N_values = [1e2, 1e3, 1e4, 1e5];

% Define mu
dx = 0.001;
mu = 4/dx^2*(sin(pi*dx/2))^2;

for N=N_values
    
    % Numerical solution
    dt = T/N;
    u_numerical = (1-dt*mu)^N;
    
    % Compute the absolute error
    error = abs(u_numerical - u_analytical)/u_analytical;
    
    % Print the error
    fprintf('%-8d %-7g %-9.3g %-5.3g\n', N, dt, error, error/dt);
end