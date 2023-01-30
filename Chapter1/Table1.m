% Compute the error of the numerical solution for solving the equation
% y'(t) = y(t), y(0) = 1, at T = 1.
clear all

% Analytical solution
T = 1;
y_analytical = exp(T);

% Print start of table
fprintf('N     dt      E         E/dt\n')

% Consider a number of different values of dt = T/N = 1/N
N_values = [5, 10, 100, 1000];

for N=N_values
    
    % Numerical solution
    dt = T/N;
    y_numerical = (1+dt)^N;
    
    % Compute the absolute error
    error = abs(y_analytical - y_numerical);
    
    % Print the error
    fprintf('%-5d %-7g %-9.3g %-5.3g\n', N, dt, error, error/dt);
end