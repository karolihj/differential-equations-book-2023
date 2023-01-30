% Run a 3D simulation of the EMI model with membrane dynamics modeled by
% the Hodgkin-Huxley model
clear all

% Set up simulation time
G.Tstop = 3;              % Total simulation time (in ms)
G.dt = 0.02;              % Time step (in ms)
Nt = round(G.Tstop/G.dt); % Number of time steps

% Set up domain
G.Lx = 0.206;     % cm 
G.Ly = 0.005;     % cm
G.Lz = 0.005;     % cm
G.dx = 0.001;     % cm
G.dy = 2.5e-4;    % cm 
G.dz = 2.5e-4;    % cm
G.Nx = round(G.Lx/G.dx)+1; % Number of grid points in the x-direction
G.Ny = round(G.Ly/G.dy)+1; % Number of grid points in the y-direction
G.Nz = round(G.Lz/G.dz)+1; % Number of grid points in the z-direction
G.N = G.Nx*G.Ny*G.Nz;      % Total number of grid points

% Set up cell geometry
G.cell_length_x = 0.2;     % Cell length (in cm)
G.cell_length_y = 0.001;   % Cell width (in cm)
G.cell_length_z = 0.001;   % Cell width (in cm)
G.cell_start_x = (G.Lx-G.cell_length_x)/2; 
G.cell_start_y = (G.Ly-G.cell_length_y)/2;
G.cell_start_z = (G.Lz-G.cell_length_z)/2;
G.stim_length = 0.05;

% Set up parameters
G.sigma_i = 4;   % mS/cm
G.sigma_e = 3;   % mS/cm
G.Cm = 1;        % uF/cm^2
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
I_ion = @(v, m, h, r) I_Na(v, m, h) + I_K(v, r) + I_L(v);

% Define rate constants
alpha_m = @(v) (gamma1*(v+gamma2))./(1-exp(-(v+gamma2)/gamma3));
beta_m = @(v) gamma4*exp(-(v+gamma5)/gamma6);
alpha_h = @(v) gamma7*exp(-(v+gamma8)/gamma9);
beta_h = @(v) gamma10./(1+exp(-(v+gamma11)/gamma12));
alpha_r = @(v) (gamma13*(v+gamma14))./(1-exp(-(v+gamma14)/gamma15));
beta_r = @(v) gamma16*exp(-(v+gamma17)/gamma18);

% Set up mesh
mesh = set_up_mesh(G);
nv = length(mesh.v);

% Set up zero rhs vector vector
d = zeros(G.N, 1);

% Set up matrix
A = set_up_matrix(G, mesh);

% Set up solution arrays
v = zeros(nv, Nt+1);
m = zeros(nv, Nt+1);
h = zeros(nv, Nt+1);
r = zeros(nv, Nt+1);
u = zeros(G.N, Nt+1);

% Define initial conditions
v(:,1) = -65;
v(mesh.to_stim,1) = -50;
m(:,1) = 0.1;
h(:,1) = 0.6;
r(:,1) = 0.3;


% Run simulation
for n = 1:Nt

    % Step 1: Solve ODE system
    v(:,n+1) = v(:,n) - G.dt/G.Cm*(I_ion(v(:,n),m(:,n),h(:,n),r(:,n)));
    m(:,n+1) = m(:,n) + G.dt*(alpha_m(v(:,n)).*(1-m(:,n)) - beta_m(v(:,n)).*m(:,n));
    h(:,n+1) = h(:,n) + G.dt*(alpha_h(v(:,n)).*(1-h(:,n)) - beta_h(v(:,n)).*h(:,n));
    r(:,n+1) = r(:,n) + G.dt*(alpha_r(v(:,n)).*(1-r(:,n)) - beta_r(v(:,n)).*r(:,n));
    
    % Step 2: Solve PDE system
    b = [d; G.Cm/G.dt*v(:,n+1)];
    X = A\b;
    u(:, n+1) = X(1:G.N);
    v(:, n+1) = X(G.N+1:end);

    % Print progress
    fprintf('%.1f%% done\n', 100*n/Nt)
end


% Select time points to plot
t_plot = [1, 2, 3];
t_idx = round(t_plot/G.dt) + 1;

% Set up figure and plot the solution
figure('Units','centimeters', 'Position', [10 10 30 22], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [27.5, 21])

z_idx = round((G.cell_start_z+G.cell_length_z)/G.dz)+1;
[x, y] = meshgrid(0:G.dx:G.Lx, 0:G.dy:G.Ly);
x=x*10; y = y*10; % Convert from cm to mm

for n=1:length(t_idx)

    % Plot the extracellular potential
    subplot(length(t_idx), 2, (n-1)*2 + 1)
    u_tmp = u(:, t_idx(n));
    u_tmp(mesh.v) = nan;
    u_tmp = reshape(u_tmp, G.Nx, G.Ny, G.Nz);
    p = pcolor(x, y, u_tmp(:,:,z_idx)');
    set(p, 'edgecolor', 'none', 'facecolor', 'interp')
    set(gca, 'fontsize', 11)
    if n==1
        title('Extracellular potential', 'fontsize', 16)
    end
    pos = get(gca, 'Position');
    pos(2) = pos(2) + 0.02*n;
    pos(3) = pos(3)+0.04;
    if n==length(t_idx)
        xlabel('x (mm)')
        c = colorbar;
        ylabel(c, 'u_e (mV)', 'fontsize', 11)
        set(c, 'Location', 'southoutside', 'fontsize', 11)
    else
        set(gca, 'XTick', [])
    end
    ylabel('y (mm)')
    clim([-0.062, 0.04])
    set(gca, 'Position', pos)
    text(-0.18, 0.5, sprintf('T = %d ms', t_plot(n)), 'units', 'normalized', ...
        'FontSize', 16, 'FontWeight', 'bold', 'Rotation', 90, ...
        'VerticalAlignment','middle', 'HorizontalAlignment','center')
    

    % Plot the membrane potential
    subplot(length(t_idx), 2, (n-1)*2 + 2)
    v_tmp = nan*ones(G.N, 1);
    v_tmp(mesh.v) = v(:, t_idx(n));
    v_tmp = reshape(v_tmp, G.Nx, G.Ny, G.Nz);
    p = pcolor(x, y, v_tmp(:,:,z_idx)');
    set(p, 'edgecolor', 'none')
    set(gca, 'fontsize', 11)
    if n==1
        title('Membrane potential', 'fontsize', 16)
    end
    pos = get(gca, 'Position');
    pos(2) = pos(2) + 0.02*n;
    pos(3) = pos(3)+0.04;
    if n==length(t_idx)
        xlabel('x (mm)')
        c = colorbar;
        ylabel(c, 'v (mV)', 'fontsize', 11)
        set(c, 'Location', 'southoutside', 'fontsize', 11)
    else
        set(gca, 'XTick', [])
    end
    set(gca, 'YTick', [])
    clim([-65, 40])
    set(gca, 'Position', pos)

end

% Save the figure
print('-dpdf', '../Figures/Ch11_Fig2.pdf')

