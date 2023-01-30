% Solve the equations of the parsimonious ventricular rabbit model using 
% an expicit numerical scheme and adjust the conductance parameters
clear all

% Set up figure
figure('Units','centimeters', 'Position', [10 10 20 12], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [18, 12])

% Set up parameters
Cm = 1;        % uF/cm^2
v_Na = 65;     % mV
v_K = -83;     % mV
b = 0.047;     % 1/mV
Em = -41;      km = -4;    
Eh = -74.9;    kh = 4.4;    tau_h_0 = 6.8;    delta_h = 0.8;

% Stimulation parameters
t_stim = 50;  % ms
d_stim = 2;   % ms
a_stim = -25; % uA/cm^2 

% Define rate constants
m_inf = @(v) 1./(1+exp((v-Em)/km));
tau_m = @(v) 0.12;
h_inf = @(v) 1./(1+exp((v-Eh)/kh));
tau_h = @(v) 2*tau_h_0*exp(delta_h*(v-Eh)/kh)./(1+exp((v-Eh)/kh));

% Set up discrerization
T = 600;          % Total simulation time (in ms)
dt = 0.001;       % Time step (in ms)
N = round(T/dt);  % Number of time steps
t = (0:dt:T);     % Time vector

adjustment_factors = [0.5, 0.75, 1, 1.25, 1.5];

for i=1:2
    for j=1:length(adjustment_factors)
        if i==1
            parameter_name = 'g_{Na}';
            g_Na = 11*adjustment_factors(j); 
            g_K = 0.3;
        elseif i==2
            parameter_name = 'g_K';
            g_Na = 11; 
            g_K = 0.3*adjustment_factors(j);
        end

        % Update legends
        legends{j} = sprintf('%g\\times%s', adjustment_factors(j), ...
            parameter_name);
        
        % Define currents
        I_Na = @(v, m, h) g_Na.*m.^3.*h.*(v-v_Na);
        I_K = @(v) g_K.*exp(-b*(v-v_K)).*(v-v_K);
        I_stim = @(t) a_stim*(t >= t_stim).*(t <= t_stim + d_stim);
        
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
        
        % Plot the upstroke
        subplot(2,3,(i-1)*3+1)
        plot(t, v, 'linewidth', 2)
        hold on
        set(gca, 'fontsize', 12)
        ylabel('v (mV)')
        xlim([51.5, 53])
        ylim([-52, 50])
        if i==1
            title('Upstroke')
        else
            xlabel('t (ms)')
        end
        
        % Plot the action potential
        subplot(2,3,(i-1)*3+2)
        plot(t, v, 'linewidth', 2)
        hold on
        set(gca, 'fontsize', 12)
        if i==1
            title('Action potential')
        else
            xlabel('t (ms)')
        end
        xlim([0, 550])
        
        % Set up legends
        subplot(2,3,(i-1)*3+3)
        plot(1, 1, 'linewidth', 2)
        hold on
        set(gca, 'fontsize', 12)
        axis off
        if j==length(adjustment_factors)
            legend(legends, 'location', 'westoutside', 'fontsize', 12)
        end
    end
end

% Save figure
print('-dpdf', '../Figures/Ch8_Fig5.pdf')