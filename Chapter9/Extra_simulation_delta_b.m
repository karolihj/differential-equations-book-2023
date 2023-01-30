% Solve the equations of the cable equation with membrane dynamics modeled
% by the Hodgkin-Huxley model using an operator splitting scheme. Compare
% the case of a constant delta to the case of a special delta at the
% boundary
clear all

% Set up parameters
Cm = 1;          % uF/cm^2
w = 0.001;       % cm
sigma_i = 4;     % mS/cm
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
I_ion = @(v, m, h, r) I_Na(v, m, h) + I_K(v, r) + I_L(v);

% Define rate constants
alpha_m = @(v) (gamma1*(v+gamma2))./(1-exp(-(v+gamma2)/gamma3));
beta_m = @(v) gamma4*exp(-(v+gamma5)/gamma6);
alpha_h = @(v) gamma7*exp(-(v+gamma8)/gamma9);
beta_h = @(v) gamma10./(1+exp(-(v+gamma11)/gamma12));
alpha_r = @(v) (gamma13*(v+gamma14))./(1-exp(-(v+gamma14)/gamma15));
beta_r = @(v) gamma16*exp(-(v+gamma17)/gamma18);

% Set up discretization
T = 7;               % Total simulation time (in ms)
dt = 0.01;           % Time step (in ms)
N = round(T/dt);     % Number of time steps
t = (0:dt:T);        % Time vector
L = 0.5;             % Length of cell/axon (in cm)
dx = 0.001;          % Spatial step (in cm)
M = round(L/dx) + 1; % Number of spatial points

% Define delta
delta = w*sigma_i/4;
delta_b = sigma_i*w*w/(w*w/dx + 4*w);

% Define the matrix
A = delta/(dx^2)*(spdiags([-1; -2*ones(M-2, 1); -1], 0, M, M) ...
    + spdiags(ones(M,1), 1, M, M) ...
    + spdiags(ones(M,1), -1, M, M));
Ab = 1/(dx^2)*(spdiags([-delta_b; -2*delta*ones(M-2, 1); -delta_b], 0, M, M) ...
    + spdiags([0; delta_b; delta*ones(M-2,1)], 1, M, M) ...
    + spdiags([delta*ones(M-2,1); delta_b; 0], -1, M, M));
I = eye(M);

%% Default case
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

% Operator splitting scheme
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
end

% Select time points to plot
t_plot = [0, 1, 3, 5, 7];
t_idx = round(t_plot/dt) + 1;
x = (0:dx:L);

% Set up figure and plot the solution
figure('Units','centimeters', 'Position', [10 10 15 24], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [15, 22])

for i=1:length(t_plot)
    subplot(length(t_plot), 1, i)
    plot(x, v(:, t_idx(i)), 'linewidth', 2)
    hold on
    set(gca, 'fontsize', 12)
    title(sprintf('t = %g ms', t_plot(i)))
    ylabel('v')
    if i == length(t_plot)
        xlabel('x (cm)')
    end
    ylim([-100, 50])
end

%% delta_b case
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

% Operator splitting scheme
for n=1:N
    % Step 1
    v(:,n+1) = (I-(dt/Cm)*Ab)\v(:, n);
    m(:,n+1) = m(:,n);
    h(:,n+1) = h(:,n);
    r(:,n+1) = r(:,n);

    % Step 2
    v(:,n+1) = v(:,n+1) - dt/Cm*(I_ion(v(:,n+1),m(:,n+1),h(:,n+1),r(:,n+1)));
    m(:,n+1) = m(:,n+1) + dt*(alpha_m(v(:,n+1)).*(1-m(:,n+1)) - beta_m(v(:,n+1)).*m(:,n+1));
    h(:,n+1) = h(:,n+1) + dt*(alpha_h(v(:,n+1)).*(1-h(:,n+1)) - beta_h(v(:,n+1)).*h(:,n+1));
    r(:,n+1) = r(:,n+1) + dt*(alpha_r(v(:,n+1)).*(1-r(:,n+1)) - beta_r(v(:,n+1)).*r(:,n+1));
end

% Select time points to plot
t_idx = round(t_plot/dt) + 1;

for i=1:length(t_plot)
    subplot(length(t_plot), 1, i)
    plot(x, v(:, t_idx(i)), ':', 'linewidth', 2)

    if i==1
        legend('\delta everywhere', '\delta_b at boundary')
    end
end

print('-dpdf', '../Figures/Ch9_extra_fig.pdf')