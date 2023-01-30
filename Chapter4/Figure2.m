% Solve the diffusion equation using an explicit and an implicit numerical 
% scheme with large time steps
clear all

% Discretization parameters
T = 0.5;
M = 51;
dx = 1/(M-1);
N = 500;
dt = T/N;
rho = dt/dx^2;

% Set up x
x = (0:dx:1);

% Select time points to plot
t_plot = [0, 0.01, 0.05, 0.1, 0.5];
t_idx = round(t_plot/dt) + 1;

% Set up figure
figure('Units','centimeters', 'Position', [10 10 25 24], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [22, 22])



%% EXPLICIT SCHEME

% Solution vector
u = zeros(M, N+1);

% Initial conditions
u(1:round(M/2), 1) = 1;
u(round(M/2)+1:end, 1) = 0;

% Numerical scheme
for n=1:N
    for j=1:M
        if j==1 % Left boundary
            u(j, n+1) = (1-2*rho)*u(j,n) + 2*rho*u(j+1,n);
        elseif j==M % Right boundary            
            u(j, n+1) = 2*rho*u(j-1,n) + (1-2*rho)*u(j,n);
        else % Inner points
            u(j, n+1) = rho*u(j-1,n) + (1-2*rho)*u(j,n) + rho*u(j+1,n);
        end
    end
end

% Plot the solution
for i=1:length(t_plot)
    subplot(length(t_plot), 2, (i-1)*2 + 1)
    plot(x, u(:, t_idx(i)), 'linewidth', 2)
    set(gca, 'fontsize', 14)
    text(0.98, 0.88, sprintf('t = %g', t_plot(i)), 'fontsize', 14, ...
        'Units', 'normalized', 'HorizontalAlignment', 'right')
    if i==1
        title('Explicit scheme', 'fontsize', 20)
    end
    ylabel('u')
    if i == length(t_plot)
        xlabel('x')
    end
end



%% IMPLICIT SCHEME

% Solution vector
u = zeros(M, N+1);

% Initial conditions
u(1:round(M/2), 1) = 1;
u(round(M/2)+1:end, 1) = 0;

% Set up matrix
A = spdiags((1+2*rho)*ones(M, 1), 0, M, M) ...
    + spdiags([0; -2*rho; -rho*ones(M-2, 1)], 1, M, M) ...
    + spdiags([-rho*ones(M-2, 1); -2*rho; 0], -1, M, M);

% Numerical scheme
for n=1:N
    u(:,n+1) = A\u(:,n);
end

% Plot the solution
for i=1:length(t_plot)
    subplot(length(t_plot), 2, (i-1)*2 + 2)
    plot(x, u(:, t_idx(i)), 'linewidth', 2)
    set(gca, 'fontsize', 14)
    text(0.98, 0.88, sprintf('t = %g', t_plot(i)), 'fontsize', 14, ...
        'Units', 'normalized', 'HorizontalAlignment', 'right')
    if i==1
        title('Implicit scheme', 'fontsize', 20)
    end
    ylabel('u')
    if i == length(t_plot)
        xlabel('x')
    end
    ylim([0,1])
end

% Save figure
print('-dpdf', '../Figures/Ch4_Fig2.pdf')