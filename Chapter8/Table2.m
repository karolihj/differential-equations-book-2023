% Solve the equations of the parsimonious ventricular rabbit model using 
% an expicit numerical scheme and print the error
clear all

% Set up parameters
Cm = 1;          % uF/cm^2
g_Na = 11;       % mS/cm^2
g_K = 0.3;       % mS/cm^2
v_Na = 65;       % mV
v_K = -83;       % mV
b = 0.047;       % 1/mV
Em = -41;      km = -4;    
Eh = -74.9;    kh = 4.4;    tau_h_0 = 6.8;    delta_h = 0.8;

% Stimulation parameters
stim_start = 0;
stim_duration = 2;

% Define currents
I_Na = @(v, m, h) g_Na.*m.^3.*h.*(v-v_Na);
I_K = @(v, r) g_K.*exp(-b*(v-v_K)).*(v-v_K);
I_stim = @(t) -25*(t >= stim_start).*(t <= stim_start + stim_duration);

% Define rate constants
m_inf = @(v) 1./(1+exp((v-Em)/km));
tau_m = @(v) 0.12;
h_inf = @(v) 1./(1+exp((v-Eh)/kh));
tau_h = @(v) 2*tau_h_0*exp(delta_h*(v-Eh)/kh)./(1+exp((v-Eh)/kh));

%% VERY FINE RESOLUTION
% Set up discrerization
T = 10;           % Total simulation time (in ms)
dt = 1e-5;        % Time step (in ms)
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

v_fine = v;

%% COARSER TIME RESOLUTIONS
% Print start of table
fprintf('dt      E         E/dt\n')

dt_values = [0.01, 0.005, 0.002, 0.001, 0.0005];

for i=1:length(dt_values)
    dt = dt_values(i);
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
    
    % Compute error
    error = abs(v(end)-v_fine(end));
    
    % Print the error
    fprintf('%-7g %-9.3g %-9.2g\n', dt, error, error/dt);

end
