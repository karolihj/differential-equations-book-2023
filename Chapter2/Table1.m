% Compute the numerical solution of the FitzHugh-Nagumo model for some
% different values for the time step dt
clear all

% Define model constants
a = -0.12;
c1 = 0.175;
c2 = 0.03;
b = 0.011;
d = 0.55;



%% Solve the system with a very fine resolution 

% Define time step
T = 5000;      % Simulation time
N = 5e6;       % Number of time points
dt = T/N;      % Time step

% Set up arrays for saving the solutions
v = zeros(N+1, 1);
w = zeros(N+1, 1);

% Define initial conditions
v(1) = 0.26; 
w(1) = 0;

% Compute the numerical solution
for n=1:N
    v(n+1) = v(n) + dt*(c1*v(n)*(v(n)-a)*(1-v(n)) - c2*w(n));
    w(n+1) = w(n) + dt*(b*(v(n)-d*w(n)));
end

dt_fine = dt;
t_fine = t;
v_fine = v;
w_fine = w;



%% Solve the system with coarser resolutions 

N_values = [500, 1000, 5000, 10000, 50000];

% Print start of table
fprintf('N       dt    E         E/dt\n')

for i=1:length(N_values)

    % Define time step
    N = N_values(i);
    dt = T/N;  

    % Set up arrays for saving the solutions
    v = zeros(N+1, 1);
    w = zeros(N+1, 1);

    % Define initial conditions
    v(1) = 0.26; 
    w(1) = 0;

    % Compute the numerical solution
    for n=1:N
        v(n+1) = v(n) + dt*(c1*v(n)*(v(n)-a)*(1-v(n)) - c2*w(n));
        w(n+1) = w(n) + dt*(b*(v(n)-d*w(n)));
    end
    
    % Compute the error
    E = abs(v_fine(end)-v(end)) + abs(w_fine(end)-w(end));
    
    % Print the error
    fprintf('%-7d %-5g %-9.3g %-8.2g\n', N, dt, E, E/dt);

end

 
