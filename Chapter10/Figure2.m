% Compute the solution of the bidomain model with membrane dynamics modeled
% by the parsimonious ventricular rabbit model
clear all

% Set up parameters
Cm = 1;           % uF/cm^2
chi = 2000;       % 1/cm
sigma_i_x = 3.0;  % mS/cm
sigma_i_y = 3.0;  % mS/cm
sigma_e_x = 10.0; % mS/cm
sigma_e_y = 10.0; % mS/cm
g_Na = 11;        % mS/cm^2
g_K = 0.3;        % mS/cm^2
v_Na = 65;        % mV
v_K = -83;        % mV
b_K = 0.047;      % 1/mV
Em = -41;      km = -4;    
Eh = -74.9;    kh = 4.4;    tau_h_0 = 6.8;    delta_h = 0.8;

% Stimulation parameters
t_stim = 0;     % Stimulation start time (in ms)
d_stim = 2;     % Stimulation current duration (in ms)
a_stim = -25;   % Stimulation current amplitude (in uA/cm^2)
l_stim = 0.25;  % Radius of the stimulation area (in cm)

% Define currents
I_Na = @(v, m, h) g_Na.*m.^3.*h.*(v-v_Na);
I_K = @(v) g_K.*exp(-b_K*(v-v_K)).*(v-v_K);
I_stim = @(t,x,y) a_stim*(t >= t_stim).*(t <= t_stim + d_stim).*(sqrt(x.^2+y.^2) <= l_stim);

% Define rate constants
m_inf = @(v) 1./(1+exp((v-Em)/km));
tau_m = @(v) 0.12;
h_inf = @(v) 1./(1+exp((v-Eh)/kh));
tau_h = @(v) 2*tau_h_0*exp(delta_h*(v-Eh)/kh)./(1+exp((v-Eh)/kh));

% Set up discrerization
T = 20;                % Total simulation time (in ms)
dt = 0.01;             % Time step (in ms)
N = round(T/dt);       % Number of time steps
t = (0:dt:T);          % Time vector
Lx = 1;                % Length of the domain in the x-direction (in cm)
Ly = 1;                % Length of the domain in the y-direction (in cm)
dx = 0.025;            % Discretization step in the x-direction (in cm)
dy = 0.025;            % Discretization step in the y-direction (in cm)
Mx = round(Lx/dx) + 1; % Number of point in the x-direction
My = round(Ly/dy) + 1; % Number of point in the y-direction
M = Mx*My;             % Total number of spatial points

% Define x and y for each global index
x = dx*rem((1:M)-1, Mx)';
y = dy*floor(((1:M)-1)/Mx)';

% Set up solution arrays
v = zeros(M, N+1);
u = zeros(M, N+1);
m = zeros(M, N+1);
h = zeros(M, N+1);

% Define initial conditions
v(:,1) = -83;
u(:,1) = 0;
m(:,1) = 0;
h(:,1) = 0.9;

% Set up matrices
Ai = set_up_matrix(Mx, My, dx, dy, sigma_i_x, sigma_i_y);
Ae = set_up_matrix(Mx, My, dx, dy, sigma_e_x, sigma_e_y);
I = eye(M);
A = [I - (dt/(chi*Cm))*Ai, -(dt/(chi*Cm))*Ai; ...
     Ai, Ai+Ae];
b = zeros(2*M, 1);

% Adjust the matrix so that ue=0 at the boundary
boundary_points = [1:Mx, (My-1)*Mx+(1:Mx), (1:My-2)*Mx+1, (2:My-1)*Mx];
for i=1:length(boundary_points)
    zero_one_vector = zeros(2*M, 1);
    zero_one_vector(M+boundary_points(i)) = 1;
    A(M+boundary_points(i),:) = zero_one_vector;
end
A = sparse(A);

% Operator splitting scheme
for n=1:N

    % Step 1: Explicit ODE
    v(:,n+1) = v(:,n) - (dt/Cm)*(I_Na(v(:,n),m(:,n),h(:,n)) + I_K(v(:,n)) + I_stim(t(n),x,y));
    m(:,n+1) = m(:,n) + dt*((m_inf(v(:,n))-m(:,n))./tau_m(v(:,n)));
    h(:,n+1) = h(:,n) + dt*((h_inf(v(:,n))-h(:,n))./tau_h(v(:,n)));

    % Step 2: Implicit PDE
    b(1:M) = v(:,n+1);      % Update right-hand side
    vu = A\b;               % Solve linear system
    v(:,n+1) = vu(1:M);     % Extract v
    u(:,n+1) = vu(M+1:2*M); % Extract u

    % Print progress
    if rem(n, round(T/(dt*100))) == 0
        fprintf('%g%% done\n', 100*n/N)
    end

end

% Select time points to plot
t_plot = [5, 10, 15, 20];
t_idx = round(t_plot/dt) + 1;
x_vec = (0:dx:Lx);
y_vec = (0:dy:Ly);
[x, y] = meshgrid(x_vec, y_vec);

% Set up figure
figure('Units','centimeters', 'Position', [10 10 22 10], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [20.5, 10])

for i=1:length(t_idx)

    % Plot the membrane potential
    subplot(2, length(t_idx), i)
    pc = pcolor(x, y, reshape(v(:,t_idx(i)), Mx, My)');
    set(pc, 'edgecolor', 'none')
    axis equal
    set(gca, 'fontsize', 11)
    caxis([-83, 30])
    title(sprintf('t = %.0f ms', t_plot(i)))

    if i==length(t_idx)
        cb = colorbar;
        ylabel(cb, 'v (mV)')
    end
    if i==1
        ylabel('y (cm)')
    end

    % Adjust subplot position
    pos = get(gca, 'Position');
    pos(1) = pos(1) - 0.065;
    pos(3) = 0.22;
    pos(4) = 0.35;
    set(gca, 'Position', pos);

    % Plot the extracellular potential
    subplot(2, length(t_idx), length(t_idx) + i)
    pc = pcolor(x, y, reshape(u(:,t_idx(i)), Mx, My)');
    set(pc, 'edgecolor', 'none')
    axis equal
    set(gca, 'fontsize', 11)
    xlabel('x (cm)')
    caxis([-20, 15])
    if i==length(t_idx)
        cb = colorbar;
        ylabel(cb, 'u_e (mV)')
    end
    if i==1
        ylabel('y (cm)')
    end

    % Adjust subplot position
    pos = get(gca, 'Position');
    pos(1) = pos(1) - 0.065;
    pos(2) = pos(2) + 0.01;
    pos(3) = 0.22;
    pos(4) = 0.35;
    set(gca, 'Position', pos);
    
end

% Save the figure
print('-dpdf', '../Figures/Ch10_Fig2.pdf')

% Compute the conduction velocity
v_th = -20;
x_start = 0.4;
x_end = 0.8;
y_start = 0.4;
y_end = 0.8;
CV = compute_cv(v, v_th, t, x_vec, y_vec, x_start, x_end, y_start, y_end);
fprintf('CV = %.1f cm/s\n', CV)


