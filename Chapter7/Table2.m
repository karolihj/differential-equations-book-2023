% Compute the integral of x^2 from 0 to 1 using two different methods
clear all

% Print start of table
fprintf('dt      E_R       E_R/dx    E_T       E_T/dx^2\n')

% Define g
g = @(x) x.^2;

% Analytical integral
L = 1;      % Length of spatial domain
I_A = 1/3;  % Analytical integral value

% Discretization values
M_values = [11, 51, 101, 501, 1001];

% Run through the M values and compute approximations of the integral
for m=1:length(M_values)
    M = M_values(m);  % Number of spatial grid points
    dx = L/(M-1);     % Distance between spatial grid points
    x = (0:dx:L);     % x values
    
    % Riemann sum
    I_R = dx*sum(g(x(1:M-1)));

    % Trapezoidal method
    I_T = dx*sum([(1/2)*g(x(1)), g(x(2:M-1)), (1/2)*g(x(M))]);

    % Compute the error
    E_R = abs(I_R-I_A);
    E_L = abs(I_T-I_A);

    % Print the error
    fprintf('%-7g %-9.3g %-9.2g %-9.3g %-9.2g\n', dx, E_R, E_R/dx, ...
        E_L, E_L/dx^2);

end