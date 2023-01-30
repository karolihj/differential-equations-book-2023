% Plot some properties related to the open probability of the parsimonious 
% ventricular rabbit model 
clear all

% Set up parameters
v_K = -83;     % mV
b = 0.047;     % 1/mV
Em = -41;      km = -4;    
Eh = -74.9;    kh = 4.4;    tau_h_0 = 6.8;    delta_h = 0.8;

% Define rate constants
m_inf = @(v) 1./(1+exp((v-Em)/km));
tau_m = @(v) 0.12;
h_inf = @(v) 1./(1+exp((v-Eh)/kh));
tau_h = @(v) 2*tau_h_0*exp(delta_h*(v-Eh)/kh)./(1+exp((v-Eh)/kh));
o_K = @(v) exp(-b*(v-v_K));

% Set up v
v = linspace(-100, 50, 201);

% Set up figure
figure('Units','centimeters', 'Position', [10 10 20 6], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [17.2, 6])

% Plot m_inf
subplot(1,3,1)
plot(v, m_inf(v), 'linewidth', 2)
set(gca, 'fontsize', 12)
title('m_{\infty}')
xlabel('v (mV)')

% Plot h_inf
subplot(1,3,2)
plot(v, h_inf(v), 'linewidth', 2)
set(gca, 'fontsize', 12)
title('h_{\infty}')
xlabel('v (mV)')

% Plot K-channel open probability
subplot(1,3,3)
plot(v, o_K(v), 'linewidth', 2)
set(gca, 'fontsize', 12)
title('e^{-b*(v-v_K)}')
xlabel('v (mV)')

% Save figure
print('-dpdf', '../Figures/Ch8_Fig3.pdf')


