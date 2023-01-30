% Compute the error of the numerical solution for an explicit, an implicit,
% and a midpoint scheme for solving the equaation 
% y'(t) = y(t), y(0) = 1, at T = 1.
clear all

% Analytical solution
T = 1;
y_analytical = exp(T);

% Print start of table
fprintf('dt      E_E       E_E/dt    E_I       E_I/dt    E_M       E_M/dt^2\n')

% Consider a number of different values of dt = T/N = 1/N
N_values = [10, 100, 1000, 10000];

for N=N_values

    % Define dt
    dt = T/N;

    % Numerical solution (explicit scheme)
    y_explicit = (1+dt)^N;

    % Numerical solution (implicit scheme)
    y_implicit = (1/(1-dt))^N;

    % Numerical solution (midpoint scheme)
    y_midpoint = ((1+dt/2)/(1-dt/2))^N;
    
    % Compute the absolute error
    error_e = abs(y_analytical - y_explicit);
    error_i = abs(y_analytical - y_implicit);
    error_m = abs(y_analytical - y_midpoint);
    
    % Print the error
    fprintf('%-7g %-9.3g %-9.3g %-9.3g %-9.3g %-9.3g %-9.3g\n', dt, error_e, error_e/dt, ...
        error_i, error_i/dt, error_m, error_m/dt^2);
end