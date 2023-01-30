% Solve the equations of the cable equation with membrane dynamics modeled
% by the Hodgkin-Huxley model using an expicit numerical scheme
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
T = 1;               % Total simulation time (in ms)
dt = 0.001;          % Time step (in ms)
N = round(T/dt);     % Number of time steps
t = (0:dt:T);        % Time vector
L = 0.5;             % Length of cell/axon (in cm)
dx = 0.001;          % Spatial step (in cm)
M = round(L/dx) + 1; % Number of spatial points

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

% Explicit numerical scheme
for n=1:N
    v(:,n+1) = (I+(dt/Cm)*A)*v(:, n) - dt/Cm*(I_ion(v(:,n),m(:,n),h(:,n),r(:,n)));
    m(:,n+1) = m(:,n) + dt*(alpha_m(v(:,n)).*(1-m(:,n)) - beta_m(v(:,n)).*m(:,n));
    h(:,n+1) = h(:,n) + dt*(alpha_h(v(:,n)).*(1-h(:,n)) - beta_h(v(:,n)).*h(:,n));
    r(:,n+1) = r(:,n) + dt*(alpha_r(v(:,n)).*(1-r(:,n)) - beta_r(v(:,n)).*r(:,n));
end

% Select time points to plot
t_plot = [0, 0.001, 0.002, 0.004, 0.005];
t_idx = round(t_plot/dt) + 1;
x = (0:dx:L);

% Set up figure and plot the solution
figure('Units','centimeters', 'Position', [10 10 12 18], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [12, 17])

for i=1:length(t_plot)
    subplot(length(t_plot), 1, i)
    plot(x, v(:, t_idx(i)), 'linewidth', 2)
    set(gca, 'fontsize', 12)
    title(sprintf('t = %g ms', t_plot(i)))
    ylabel('v (mV)')
    if i == length(t_plot)
        xlabel('x (cm)')
    end
    xlim([0, 0.1])
    ylim([min(-100, min(v(:, t_idx(i)))), max(50, max(v(:, t_idx(i))))])
end

% Save the figure
print('-dpdf', '../Figures/Ch9_Fig2.pdf')
