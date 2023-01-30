% Solve the equations of the Hodgkin-Huxley model using an expicit
% numerical scheme. Compare the solution to the solution with a very fine
% time step
clear all

% Set up parameters
Cm = 1;          % uF/cm^2
g_Na = 120;      % mS/cm^2
g_K = 36;        % mS/cm^2
g_L = 0.3;       % mS/cm^2
v_Na = 50;       % mV
v_K = -77;       % mV
v_L = -54.4;     % mV
gamma1 = 0.1;     gamma2 = 40;    gamma3 = 10;     
gamma4 = 4;       gamma5 = 65;    gamma6 = 18;    
gamma7 = 0.07;    gamma8 = 65;    gamma9 = 20;     
gamma10 = 1;      gamma11 = 35;   gamma12 = 10;   
gamma13 = 0.01;   gamma14 = 55;   gamma15 = 10;    
gamma16 = 0.125;  gamma17 = 65;   gamma18 = 80;

% Define currents
I_Na = @(v, m, h) g_Na.*m.^3.*h.*(v-v_Na);
I_K = @(v, r) g_K.*r.^4.*(v-v_K);
I_L = @(v) g_L.*(v-v_L);

% Define rate constants
alpha_m = @(v) (gamma1*(v+gamma2))./(1-exp(-(v+gamma2)/gamma3));
beta_m = @(v) gamma4*exp(-(v+gamma5)/gamma6);
alpha_h = @(v) gamma7*exp(-(v+gamma8)/gamma9);
beta_h = @(v) gamma10./(1+exp(-(v+gamma11)/gamma12));
alpha_r = @(v) (gamma13*(v+gamma14))./(1-exp(-(v+gamma14)/gamma15));
beta_r = @(v) gamma16*exp(-(v+gamma17)/gamma18);

%% SOLVE SYSTEM WITH VERY FINE RESOLUTION

% Set up discrerization
T = 3;            % Total simulation time (in ms)
dt = 1e-6;        % Time step (in ms)
N = round(T/dt);  % Number of time steps

% Set up solution vectors
v = zeros(N+1, 1);
m = zeros(N+1, 1);
h = zeros(N+1, 1);
r = zeros(N+1, 1);

% Define initial conditions
v(1) = -60;
m(1) = 0.1;
h(1) = 0.6;
r(1) = 0.3;

% Explicit numerical scheme
for n=1:N
    v(n+1) = v(n) - (dt/Cm)*(I_Na(v(n),m(n),h(n)) + I_K(v(n),r(n)) + I_L(v(n)));
    m(n+1) = m(n) + dt*(alpha_m(v(n))*(1-m(n)) - beta_m(v(n))*m(n));
    h(n+1) = h(n) + dt*(alpha_h(v(n))*(1-h(n)) - beta_h(v(n))*h(n));
    r(n+1) = r(n) + dt*(alpha_r(v(n))*(1-r(n)) - beta_r(v(n))*r(n));
end

v_fine = v;

%% SOLVE SYSTEM WITH LARGER VALUES OF DT
dt_values = [0.01, 0.005, 0.001, 0.0005, 0.0001];

% Print start of table
fprintf('dt      E         E/dt\n')

for i=1:length(dt_values)
    dt = dt_values(i);
    N = round(T/dt);  % Number of time steps

    % Set up solution vectors
    v = zeros(N+1, 1);
    m = zeros(N+1, 1);
    h = zeros(N+1, 1);
    r = zeros(N+1, 1);

    % Define initial conditions
    v(1) = -60;
    m(1) = 0.1;
    h(1) = 0.6;
    r(1) = 0.3;

    % Explicit numerical scheme
    for n=1:N
        v(n+1) = v(n) - (dt/Cm)*(I_Na(v(n),m(n),h(n)) + I_K(v(n),r(n)) + I_L(v(n)));
        m(n+1) = m(n) + dt*(alpha_m(v(n))*(1-m(n)) - beta_m(v(n))*m(n));
        h(n+1) = h(n) + dt*(alpha_h(v(n))*(1-h(n)) - beta_h(v(n))*h(n));
        r(n+1) = r(n) + dt*(alpha_r(v(n))*(1-r(n)) - beta_r(v(n))*r(n));
    end

    % Compute error
    error = abs(v_fine(end)-v(end));

    % Print the error
    fprintf('%-7g %-9.3g %-9.2g\n', dt, error, error/dt);
end
