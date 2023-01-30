% Solve the equations of the Hodgkin-Huxley model using an expicit
% numerical scheme
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

% Set up discrerization
T = 10;           % Total simulation time (in ms)
dt = 0.001;       % Time step (in ms)
N = round(T/dt);  % Number of time steps
t = (0:dt:T);     % Time vector

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

% Set up figure and plot the solution
figure('Units','centimeters', 'Position', [10 10 20 25], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [19, 22])

subplot(4,2,1)
plot(t, v, 'linewidth', 2)
set(gca, 'fontsize', 12)
title('v')
ylabel('v (mV)')
ylim([-80, 50])

subplot(4,2,2)
plot(t, m, 'linewidth', 2)
set(gca, 'fontsize', 12)
title('m')
ylabel('m')
ylim([0, 1])

subplot(4,2,3)
plot(t, h, 'linewidth', 2)
set(gca, 'fontsize', 12)
title('h')
ylabel('h')
ylim([0, 1])

subplot(4,2,4)
plot(t, r, 'linewidth', 2)
set(gca, 'fontsize', 12)
title('r')
ylabel('r')
ylim([0, 1])

subplot(4,2,5)
plot(t, I_Na(v, m, h), 'linewidth', 2)
set(gca, 'fontsize', 12)
title('I_{Na}')
ylabel('I_{Na} (\muA/cm^2)')

subplot(4,2,6)
plot(t, I_K(v, r), 'linewidth', 2)
set(gca, 'fontsize', 12)
title('I_{K}')
ylabel('I_{K} (\muA/cm^2)')
xlabel('t (ms)')

subplot(4,2,7)
plot(t, I_L(v), 'linewidth', 2)
set(gca, 'fontsize', 12)
title('I_{L}')
ylabel('I_{L} (\muA/cm^2)')
xlabel('t (ms)')


% Save figure
print('-dpdf', '../Figures/Ch8_Fig1.pdf')


