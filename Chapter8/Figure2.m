% Solve the equations of the Hodgkin-Huxley model using an expicit
% numerical scheme. Investigate how the solutions are affected by 
% adjusting g_Na, g_K, and g_L
clear all

% Set up figure
figure('Units','centimeters', 'Position', [10 10 20 15], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [18, 15])

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

for j=1:3
    for i=1:3
        if j==1 
            % Adjust g_Na
            adjustment_factors = [0.8, 1, 1.5];
            parameter_name = 'g_{Na}';
            g_Na = 120*adjustment_factors(i); 
            g_K = 36;
            g_L = 0.3;
        elseif j==2
            % Adjust g_K
            adjustment_factors = [0.75, 1, 1.25];
            parameter_name = 'g_{K}';
            g_Na = 120;
            g_K = 36*adjustment_factors(i); 
            g_L = 0.3;
        else
            % Adjust g_L
            adjustment_factors = [0.5, 1, 1.5];
            parameter_name = 'g_{L}';
            g_Na = 120;
            g_K = 36; 
            g_L = 0.3*adjustment_factors(i);
        end

        % Update legends
        legends{i} = sprintf('%g\\times%s', adjustment_factors(i), ...
            parameter_name);

        % Define currents
        I_Na = @(v, m, h) g_Na.*m.^3.*h.*(v-v_Na);
        I_K = @(v, r) g_K.*r.^4.*(v-v_K);
        I_L = @(v) g_L.*(v-v_L);
        
        % Set up solution vectors
        v = zeros(N+1, 1);
        m = zeros(N+1, 1);
        h = zeros(N+1, 1);
        r = zeros(N+1, 1);
        t = (0:dt:T);
        
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

        % Adjust time for plotting
        [~, t_idx] = max((v(2:end)-v(1:end-1))/dt);
        t_idx = t_idx- round(1/dt);
        t = t - t(t_idx);
        t = t(t_idx:end);
        v = v(t_idx:end);

        % Plot the upstroke
        subplot(3,3,(j-1)*3+1)
        plot(t, v, 'linewidth', 2)
        hold on
        set(gca, 'fontsize', 12)
        ylabel('v (mV)')
        if j==1
            title('Upstroke')
        elseif j==3
            xlabel('t (ms)')
        end
        ylim([-80, 50])
        xlim([0.5, 1.5])

        % Plot the action potential
        subplot(3,3,(j-1)*3+2)
        plot(t, v, 'linewidth', 2)
        hold on
        set(gca, 'fontsize', 12)
        if j==1
            title('Action potential')
        elseif j==3
            xlabel('t (ms)')
        end
        ylim([-80, 50])
        xlim([0, 4])

        % Set up legends
        subplot(3,3,(j-1)*3+3)
        plot(1, 1, 'linewidth', 2)
        hold on
        set(gca, 'fontsize', 12)
        axis off
        if i==length(adjustment_factors)
            legend(legends, 'location', 'westoutside', 'fontsize', 12)
        end

    end
end


% Save figure
print('-dpdf', '../Figures/Ch8_Fig2.pdf')