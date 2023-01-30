function mesh = set_up_mesh(G)
%set_up_mesh(G) Set up vectors containing the indices of the different 
%types of nodes

% Read geometry
Nx = G.Nx;
Ny = G.Ny;
N = G.N;

% Set up indices for outer boundary
mesh.e_sw = 1;
mesh.e_se = Nx;
mesh.e_nw = (Ny-1)*Nx + 1;
mesh.e_ne = N;
mesh.e_w = ((2:Ny-1)-1)*Nx + 1;
mesh.e_e = ((2:Ny-1)-1)*Nx + Nx;
mesh.e_s = 2:Nx-1;
mesh.e_n = (Ny-1)*Nx + (2:Nx-1);

% Set up indices for the remaining nodes
indices = ones(1, N);
indices([mesh.e_ne, mesh.e_sw, mesh.e_se, mesh.e_nw, mesh.e_w, ...
    mesh.e_e, mesh.e_s, mesh.e_n]) = 0;
mesh.e = round(find(indices));


end