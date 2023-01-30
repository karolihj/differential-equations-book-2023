% Compute the numerical solution of the FitzHugh-Nagumo model
clear all

% Define model constants
a = -0.12;
c1 = 0.175;
c2 = 0.03;
b = 0.011;
d = 0.55;

% Define time step
T = 5000;     % Simulation time
N = 5000;     % Number of time points
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

% Set up figure and plot the solutions
figure('Units','centimeters', 'Position', [10 10 26 9], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [24, 9])

subplot(1,2,1)
plot(t, v, 'linewidth', 2)
set(gca, 'fontsize', 14)
xlabel('t')
ylabel('v')
title('v', 'fontsize', 18)

subplot(1,2,2)
plot(t, w, 'linewidth', 2)
set(gca, 'fontsize', 14)
xlabel('t')
ylabel('w')
title('w', 'fontsize', 18)

% Save figure
print('-dpdf', '../Figures/Ch2_Fig1.pdf')



% Set up figure and plot the solutions in a parametric plot
figure('Units','centimeters', 'Position', [10 10 12 10], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [12, 10])
plot(v, w, 'linewidth', 2)
set(gca, 'fontsize', 14)
xlabel('v')
ylabel('w')

% Save figure
print('-dpdf', '../Figures/Ch2_Fig2.pdf')


