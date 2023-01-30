function CV = solve_system_and_compute_cv(g_Na, g_K, g_L, w, sigma_i, ...
    sim_number, sim_tot)
%CV = solve_system_and_compute_cv(g_Na, g_K, g_L, w, sigma_i, sim_number, sim_tot)

% Set up parameters
Cm = 1;          % uF/cm^2
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
I_ion = @(v, m, h, r) I_Na(v, m, h) + I_K(v, r) + I_L(v);

% Define rate constants
alpha_m = @(v) (gamma1*(v+gamma2))./(1-exp(-(v+gamma2)/gamma3));
beta_m = @(v) gamma4*exp(-(v+gamma5)/gamma6);
alpha_h = @(v) gamma7*exp(-(v+gamma8)/gamma9);
beta_h = @(v) gamma10./(1+exp(-(v+gamma11)/gamma12));
alpha_r = @(v) (gamma13*(v+gamma14))./(1-exp(-(v+gamma14)/gamma15));
beta_r = @(v) gamma16*exp(-(v+gamma17)/gamma18);

% Set up discretization
T = 15;              % Total simulation time (in ms)
dt = 0.02;           % Time step (in ms)
t = (0:dt:T);        % Time vector
N = round(T/dt);     % Number of time steps
L = 0.5;             % Length of cell/axon (in cm)
dx = 0.001;          % Spatial step (in cm)
M = round(L/dx) + 1; % Number of spatial points
x = (0:dx:L);        % Spatial points

% Define delta
delta = w*sigma_i/4;

% Define the matrix
A = delta/(dx^2)*(spdiags([-1; -2*ones(M-2, 1); -1], 0, M, M) ...
    + spdiags(ones(M,1), 1, M, M) ...
    + spdiags(ones(M,1), -1, M, M));
I = eye(M);


% Set up solution vectors
v = zeros(M, N+1);
m = zeros(M, N+1);
h = zeros(M, N+1);
r = zeros(M, N+1);

% Define initial conditions
v(:,1) = -65;
v(1:round(0.05/dx)+1,1) = -50;
m(:,1) = 0.1;
h(:,1) = 0.6;
r(:,1) = 0.3;

% Operator splitting scheme scheme
for n=1:N

    % Step 1
    v(:,n+1) = (I-(dt/Cm)*A)\v(:, n);
    m(:,n+1) = m(:,n);
    h(:,n+1) = h(:,n);
    r(:,n+1) = r(:,n);

    % Step 2
    v(:,n+1) = v(:,n+1) - dt/Cm*(I_ion(v(:,n+1),m(:,n+1),h(:,n+1),r(:,n+1)));
    m(:,n+1) = m(:,n+1) + dt*(alpha_m(v(:,n+1)).*(1-m(:,n+1)) - beta_m(v(:,n+1)).*m(:,n+1));
    h(:,n+1) = h(:,n+1) + dt*(alpha_h(v(:,n+1)).*(1-h(:,n+1)) - beta_h(v(:,n+1)).*h(:,n+1));
    r(:,n+1) = r(:,n+1) + dt*(alpha_r(v(:,n+1)).*(1-r(:,n+1)) - beta_r(v(:,n+1)).*r(:,n+1));

    % Print progress
    if rem(n-1, 1000) == 0
        fprintf('%.0f%% done\n', 100*(sim_number-1)/sim_tot + 100*n/N/sim_tot)
    end
end

% Compute cv
x_start = 0.2;
x_end = 0.4;
v_th = 0;
CV = compute_cv(v, v_th, t, x, x_start, x_end);

end