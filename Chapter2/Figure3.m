% Compute the numerical solution of the FitzHugh-Nagumo model for some
% different values for the time step dt and plot the solution
clear all

% Set up figure
figure('Units','centimeters', 'Position', [10 10 25 25], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [22, 23])

% Define model constants
a = -0.12;
c1 = 0.175;
c2 = 0.03;
b = 0.011;
d = 0.55;


%% Solve the system with a very fine resolution        

% Define time step
T = 5000;     % Simulation time
N = 5e6;      % Number of time points
dt = T/N;     % Time step
t = (0:dt:T); % Time vector

% Set up arrays for saving the solutions
v = zeros(N+1, 1);
w = zeros(N+1, 1);

% Define initial conditions
v(1) = 0.26; 
w(1) = 0;

% Compute the numerical solution
for n=1:N
    v(n+1) = v(n) + dt*(c1*v(n)*(v(n)-a)*(1-v(n)) - c2*w(n));
    w(n+1) = w(n) + dt*(b*(v(n)-d*w(n)));
end

dt_fine = dt;
t_fine = t;
v_fine = v;
w_fine = w;




%% Solve the system with coarser resolutions         

N_values = [500, 1000, 5000];

for i=1:length(N_values)

    % Define time step
    N = N_values(i);
    dt = T/N;         % Time step
    t = (0:dt:T);     % Time vector

    % Set up arrays for saving the solutions
    v = zeros(N+1, 1);
    w = zeros(N+1, 1);

    % Define initial conditions
    v(1) = 0.26; 
    w(1) = 0;

    % Compute the numerical solution
    for n=1:N
        v(n+1) = v(n) + dt*(c1*v(n)*(v(n)-a)*(1-v(n)) - c2*w(n));
        w(n+1) = w(n) + dt*(b*(v(n)-d*w(n)));
    end
    
    % Plot v solution
    subplot(3,2,(i-1)*2+1)
    plot(t_fine, v_fine, 'linewidth', 2)
    hold on
    plot(t, v, ':', 'linewidth', 2)
    set(gca, 'fontsize', 14)
    if i==1
        title('v', 'fontsize', 18)
    elseif i == 3
        xlabel('t')
    end
    ylabel('v')
    xlim([4000, 5000])
    ylim([-0.5, 1.5])
    legend({sprintf('\\Deltat = %g', dt_fine), ...
        sprintf('\\Deltat = %g', dt)})
   
    % Plot w solution
    subplot(3,2,(i-1)*2+2)
    plot(t_fine, w_fine, 'linewidth', 2)
    hold on
    plot(t, w, ':', 'linewidth', 2)
    set(gca, 'fontsize', 14)
    if i==1
        title('w', 'fontsize', 18)
    elseif i == 3
        xlabel('t')
    end
    ylabel('w')
    xlim([4000, 5000])
    ylim([-0.5, 1.5])
    legend({sprintf('\\Deltat = %g', dt_fine), ...
        sprintf('\\Deltat = %g', dt)})
    
end

% Save figure
print('-dpdf', '../Figures/Ch2_Fig3.pdf')
 
