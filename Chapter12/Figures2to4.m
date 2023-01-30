% Run a 2D simulation of the PNP model with a K+ channel embedded
% in the membrane separating the intracellular and extracellular spaces
clear all

% Set up geometry parameters
G.Li = 50;        % Length of intracellular domain (in nm)
G.Lm = 5;         % Membrane width (in nm)
G.Le = 50;        % Length of extracellular domain (in nm)
G.Ly = 50;        % Length of domain in the y-direction (in nm)
G.w_channel = 5;  % Width of the potassium channel (in nm)
G.dx = 0.5;       % Distance between spatial mesh points (in nm)
G.dy = 0.5;       % Distance between spatial mesh points (in nm=
G.Lx = G.Li + G.Lm + G.Le;  % Total length of the domain in the x-direction
G.channel_start = (G.Ly-G.w_channel)/2;  % y-value for channel to start
G.channel_end = G.channel_start + G.w_channel;  % y-value for channel to end
G.Nx = round(G.Lx/G.dx);  % Number of mesh points in the x-direction
G.Ny = round(G.Ly/G.dy);  % Number of mesh points in the y-direction
G.N = G.Nx*G.Ny;          % Total number of mesh points

% Set up temporal discterization parameters
G.dt = 2e-8;                % Time step (in ms)
G.Tstop = 5e-5;             % Total simulation time (in ms)
G.Nt = round(G.Tstop/G.dt); % Number of time steps

% Physical constants
G.num_c = 4;                  % Number of ion species
G.eps_r = 80;                 % Relative permittivity in O_i and O_e
G.eps_m = 2;                  % Relative permittivity in O_m and O_c
G.eps_0 = 8.854e3;            % Vacuum permittivity (in fF/m)
G.eps = G.eps_r*G.eps_0;      % Permittivity in O_i and O_e (in fF/m)
G.eps_mem = G.eps_m*G.eps_0;  % Permittivity in O_m and O_c (in fF/m)
G.eps_func = @permittivity;   % Function for the spatially varying permittivity
G.rho0 = 0;                   % Default rho0
G.e = 1.60217662e-19;         % Elementary charge (in C)
G.kB = 1.38064852e-20;        % Boltzmann constant (in mJ/K)
G.T = 310;                    % Temperature (in K) 
G.D = [1.33e6; 1.96e6; 0.71e6; 2.03e6];     % Diffusion coefficient for Na, K, Ca, Cl (in nm^2/ms)
G.D_func = @diffusion_coefficient;          % Function for the spatiallt varying diffusion coefficient
G.c0e = [100; 5; 1.4; 107.8];               % Intracellular initial concentrations, [Na], [K], [Ca], [Cl] (in mM)
G.c0i = [12; 125; 0.0001; 137.0002];        % Extracellular initial concentrations, [Na], [K], [Ca], [Cl] (in mM)
G.initial_conditions = @initial_conditions; % Function for the spatially varying initial conditions
G.z = [1, 1, 2, -1];                        % Valence for Na, K, Ca, Cl
G.F = 96485.3365;                           % Faraday's constant (in C/mol)

% Set up mesh
mesh = set_up_mesh(G);

% Set up matrix for the phi system
A = set_up_phi_matrix(G, mesh);

% Set up matrix for the concentration system that does not depend on phi
B1 = set_up_C_matrix(G, mesh, 1);
B2 = set_up_C_matrix(G, mesh, 2);
B3 = set_up_C_matrix(G, mesh, 3);
B4 = set_up_C_matrix(G, mesh, 4);

% Set up initial conditions
c = zeros(G.N, G.num_c);
[x, y] = x_y_from_idx(1:G.N, G);
for i=1:G.num_c
    c(:, i) = G.initial_conditions(x, y, G, i);
end

% Update rho0 so that the domain is electroneutral at t=0
for i=1:G.num_c
    G.rho0 = G.rho0 - G.F*G.z(i)*c(:,i);
end

% Set up matrices for saving the solution
phi = zeros(G.N, G.Nt+1);
C = zeros(G.N, G.num_c, G.Nt+1);
C(:,:,1) = c;

% Run simulation
for n = 1:G.Nt
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Step 1: Solve system for the electric potential %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Update rhs-vector
    b = G.rho0;
    for i=1:G.num_c
        b = b + G.F*G.z(i)*C(:,i,n); 
    end
    b = -b;
    b([mesh.e_ne, mesh.e_e, mesh.e_se]) = 0; % Dirichlet boundary condition
    
    % Solve system
    phi(:,n+1) = A\b;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   Step 2: Solve system for the concentrations   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Update the part of the matrix that depends on phi
    B1phi = set_up_Cc_matrix(G, mesh, 1, phi(:,n+1));
    B2phi = set_up_Cc_matrix(G, mesh, 2, phi(:,n+1));
    B3phi = set_up_Cc_matrix(G, mesh, 3, phi(:,n+1));
    B4phi = set_up_Cc_matrix(G, mesh, 4, phi(:,n+1));
    
    % Update concentrations
    C(:, 1, n+1) = (B1+B1phi)\C(:,1,n);
    C(:, 2, n+1) = (B2+B2phi)\C(:,2,n);
    C(:, 3, n+1) = (B3+B3phi)\C(:,3,n);
    C(:, 4, n+1) = (B4+B4phi)\C(:,4,n);
    
    % Print remaining time
    fprintf('%.1f%% done\n', 100*n/G.Nt)

end

%% Figure 2: Plot the membrane potential
% Set up figure
figure('Units','centimeters', 'Position', [10 10 12 8], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [12, 8])

% Set up time vector
t = (0:G.dt:G.Tstop)*1e6; % Convert time from ms to ns


% Plot v
i_idx = round(G.Li/G.dx);
e_idx = i_idx + round(G.Lm/G.dx) + 1;
plot(t, phi(i_idx,:)-phi(e_idx), 'linewidth', 2)
set(gca, 'fontsize', 12)
ylabel('v (mV)')
title('v = \phi_i - \phi_e', 'fontsize', 16)
xlabel('t (ns)')
ylim([-85, 0])


% Save figure
print('-dpdf', '../Figures/Ch12_Fig2.pdf')


%% Figure 3: Plot the solution in the intracellular and extracellular 
%% parts of the domain separately 
% Set up figure
figure('Units','centimeters', 'Position', [10 10 28 15], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [25.5, 14.2])

% Set up meshgrid for plotting
[x, y] = meshgrid((G.dx/2:G.dx:G.Lx-G.dx/2)', (G.dy/2:G.dy:G.Ly-G.dy/2)');

% Plot phi (intracellular)
subplot(2, 5, 1)
u = reshape(phi(:,end), G.Nx, G.Ny);
pc = pcolor(x, y, u');
set(pc, 'edgecolor', 'none')
xlim([G.Li-5, G.Li])
set(gca, 'fontsize', 12)
ylabel('y (nm)')
title('\phi', 'fontsize', 16)
c = colorbar;
ylabel(c, 'mV')
clim([-83.05, -82.8])
text(-0.9, 0.5, '\Omega_i', 'FontSize', 22, 'FontWeight', 'bold', ...
    'Rotation', 90, 'VerticalAlignment', 'middle', ...
    'HorizontalAlignment','center', 'units', 'normalized')

% Plot phi (extracellular)
subplot(2, 5, 5+1)
u = reshape(phi(:,end), G.Nx, G.Ny);
pc = pcolor(x, y, u');
set(pc, 'edgecolor', 'none')
xlim([G.Li+G.Lm, G.Li+G.Lm+5])
set(gca, 'fontsize', 12)
ylabel('y (nm)')
xlabel('x (nm)')
c = colorbar;
ylabel(c, 'mV')
clim([-0.4, 0])
text(-0.9, 0.5, '\Omega_e', 'FontSize', 22, 'FontWeight', 'bold', ...
    'Rotation', 90, 'VerticalAlignment', 'middle', ...
    'HorizontalAlignment','center', 'units', 'normalized')

i_lims = [11.9, 12.01; 124, 124.95; 0.000098, 0.0001002; 136.95, 137.8];
e_lims = [99.8, 101.25; 5.1, 5.35; 1.395, 1.43; 106.2, 108];
titles = {'Na^+', 'K^+', 'Ca^{2+}', 'Cl^{-}'};

for i=1:4

    % Plot concentration (intracellular)
    subplot(2,5,1+i)
    u = reshape(C(:,i,end), G.Nx, G.Ny);
    pc = pcolor(x, y, u');
    set(pc, 'edgecolor', 'none')
    set(gca, 'fontsize', 12, 'YTick', [])
    xlim([G.Li-5, G.Li])
    title(titles{i}, 'fontsize', 16)
    c = colorbar;
    ylabel(c, 'mM')
    clim(i_lims(i,:))

    % Plot concentration (extracellular)
    subplot(2,5,6+i)
    u = reshape(C(:,i,end), G.Nx, G.Ny);
    pc = pcolor(x, y, u');
    set(pc, 'edgecolor', 'none')
    xlabel('x (nm)')
    set(gca, 'fontsize', 12, 'YTick', [])
    xlim([G.Li+G.Lm, G.Li+G.Lm+5])
    c = colorbar;
    ylabel(c, 'mM')
    clim(e_lims(i,:))
end

% Save figure
print('-dpdf', '../Figures/Ch12_Fig3.pdf')

%% Figure 4: Plot the solution along lines in the intracellular and 
%% extracellular parts of the domain
% Set up figure
figure('Units','centimeters', 'Position', [10 10 38 12], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [33, 12])

x = (G.dx/2:G.dx:G.Lx-G.dx/2)';

% Plot phi (intracellular)
subplot(2, 5, 1)
u = reshape(phi(:,end), G.Nx, G.Ny);
plot(x, u(:, 1), x, u(:, round(G.Ny/2)), 'linewidth', 2);
xlim([G.Li-5, G.Li-G.dx/2])
set(gca, 'fontsize', 13)
ylabel('mV')
title('\phi', 'fontsize', 18)
legend('y = 0 nm', 'y = 25 nm')
text(-0.65, 0.5, '\Omega_i', 'FontSize', 24, 'FontWeight', 'bold', ...
    'Rotation', 90, 'VerticalAlignment', 'middle', ...
    'HorizontalAlignment','center', 'units', 'normalized')

% Plot phi (extracellular)
subplot(2, 5, 5+1)
u = reshape(phi(:,end), G.Nx, G.Ny);
plot(x, u(:, 1), x, u(:, round(G.Ny/2)), 'linewidth', 2);
xlim([G.Li+G.Lm+G.dx/2, G.Li+G.Lm+5])
set(gca, 'fontsize', 13)
ylabel('mV')
xlabel('x (nm)')
text(-0.65, 0.5, '\Omega_e', 'FontSize', 24, 'FontWeight', 'bold', ...
    'Rotation', 90, 'VerticalAlignment', 'middle', ...
    'HorizontalAlignment','center', 'units', 'normalized')


for i=1:4

    % Plot concentration (intracellular)
    subplot(2,5,1+i)
    u = reshape(C(:,i,end), G.Nx, G.Ny);
    plot(x, u(:, 1), x, u(:, round(G.Ny/2)), 'linewidth', 2);
    set(gca, 'fontsize', 13)
    xlim([G.Li-5, G.Li-G.dx/2])
    title(titles{i}, 'fontsize', 18)
    ylabel('mM')

    % Plot concentration (extracellular)
    subplot(2,5,6+i)
    u = reshape(C(:,i,end), G.Nx, G.Ny);
    plot(x, u(:, 1), x, u(:, round(G.Ny/2)), 'linewidth', 2);
    xlabel('x (nm)')
    set(gca, 'fontsize', 13)
    xlim([G.Li+G.Lm+G.dx/2, G.Li+G.Lm+5])
    ylabel('mM')
end

% Save the figure
print('-dpdf', '../Figures/Ch12_Fig4.pdf')