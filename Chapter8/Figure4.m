% Solve the equations of the parsimonious ventricular rabbit model using 
% an expicit numerical scheme
clear all

% Set up parameters
Cm = 1;        % uF/cm^2
g_Na = 11;     % mS/cm^2
g_K = 0.3;     % mS/cm^2
v_Na = 65;     % mV
v_K = -83;     % mV
b = 0.047;     % 1/mV
Em = -41;      km = -4;    
Eh = -74.9;    kh = 4.4;    tau_h_0 = 6.8;    delta_h = 0.8;

% Stimulation parameters
t_stim = 50;
d_stim = 2;
a_stim = -25;

% Define currents
I_Na = @(v, m, h) g_Na.*m.^3.*h.*(v-v_Na);
I_K = @(v) g_K.*exp(-b*(v-v_K)).*(v-v_K);
I_stim = @(t) a_stim*(t >= t_stim).*(t <= t_stim + d_stim);

% Define rate constants
m_inf = @(v) 1./(1+exp((v-Em)/km));
tau_m = @(v) 0.12;
h_inf = @(v) 1./(1+exp((v-Eh)/kh));
tau_h = @(v) 2*tau_h_0*exp(delta_h*(v-Eh)/kh)./(1+exp((v-Eh)/kh));

% Set up discrerization
T = 400;          % Total simulation time (in ms)
dt = 0.001;       % Time step (in ms)
N = round(T/dt);  % Number of time steps
t = (0:dt:T);     % Time vector

% Set up solution vectors
v = zeros(N+1, 1);
m = zeros(N+1, 1);
h = zeros(N+1, 1);

% Define initial conditions
v(1) = -83;
m(1) = 0;
h(1) = 0.9;

% Explicit numerical scheme
for n=1:N
    v(n+1) = v(n) - (dt/Cm)*(I_Na(v(n),m(n),h(n)) + I_K(v(n)) + I_stim(t(n)));
    m(n+1) = m(n) + dt*((m_inf(v(n))-m(n))/tau_m(v(n)));
    h(n+1) = h(n) + dt*((h_inf(v(n))-h(n))/tau_h(v(n)));
end

% Set up figure and plot the solution
figure('Units','centimeters', 'Position', [10 10 20 12], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [17.2, 12])
subplot(2,3,1)
plot(t, v, 'linewidth', 2)
set(gca, 'fontsize', 12)
title('v')
ylabel('v (mV)')

subplot(2,3,2)
plot(t, m, 'linewidth', 2)
set(gca, 'fontsize', 12)
title('m')
ylabel('m')
ylim([0, 1])

subplot(2,3,3)
plot(t, h, 'linewidth', 2)
set(gca, 'fontsize', 12)
title('h')
ylabel('h')
ylim([0, 1])

subplot(2,3,4)
plot(t, I_Na(v, m, h), 'linewidth', 2)
set(gca, 'fontsize', 12)
title('I_{Na}')
ylabel('I_{Na} (\muA/cm^2)')
xlabel('t (ms)')

subplot(2,3,5)
plot(t, I_Na(v, m, h), 'linewidth', 2)
set(gca, 'fontsize', 12)
title('I_{Na}')
ylabel('I_{Na} (\muA/cm^2)')
xlabel('t (ms)')
xlim([51, 53])


subplot(2,3,6)
plot(t, I_K(v), 'linewidth', 2)
set(gca, 'fontsize', 12)
title('I_{K}')
ylabel('I_{K} (\muA/cm^2)')
xlabel('t (ms)')

% Save figure
print('-dpdf', '../Figures/Ch8_Fig3.4df')


