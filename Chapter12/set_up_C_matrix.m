function A = set_up_C_matrix(G, mesh, c_idx)
%set_up_C_matrix Set up matrix for the finite difference equations for the
%concentrations (except for the part that depends on phi)

% Load parameters
N = G.N;
Nx = G.Nx;
dt = G.dt;
dx = G.dx;
dy = G.dy;
D_func = G.D_func;


% Load grid point types
e_ne = mesh.e_ne;
e_sw = mesh.e_sw;
e_se = mesh.e_se;
e_nw = mesh.e_nw;
e_w = mesh.e_w;
e_e = mesh.e_e;
e_s = mesh.e_s;
e_n = mesh.e_n;
e = mesh.e;

vec = zeros(N, 1);
vec_kp = zeros(N, 1);
vec_km = zeros(N, 1);
vec_jp = zeros(N, 1);
vec_jm = zeros(N, 1);


%%%%%%% BOUNDARY %%%%%%%
    
% Neumann boundary conditions for south and north boundary
% 1a) Set up rows for the south boundary
index = e_s;
[x, y] = x_y_from_idx(index, G);
vec(index) = 1 + dt*(D_func(x+dx/2, y, G, c_idx)/(dx*dx) + ...
    D_func(x-dx/2, y, G, c_idx)/(dx*dx) + ...
    D_func(x, y+dy/2, G, c_idx)/(dy*dy));
vec_kp(index+1) = -dt*D_func(x+dx/2, y, G, c_idx)/(dx*dx);
vec_km(index-1) = -dt*D_func(x-dx/2, y, G, c_idx)/(dx*dx);
vec_jp(index+Nx) = -dt*D_func(x, y+dy/2, G, c_idx)/(dy*dy);


% 2a) Set up rows for the north boundary
index = e_n;
[x, y] = x_y_from_idx(index, G);
vec(index) = 1 + dt*(D_func(x+dx/2, y, G, c_idx)/(dx*dx) + ...
    D_func(x-dx/2, y, G, c_idx)/(dx*dx) + ...
    D_func(x, y-dy/2, G, c_idx)/(dy*dy));
vec_kp(index+1) = -dt*D_func(x+dx/2, y, G, c_idx)/(dx*dx);
vec_km(index-1) = -dt*D_func(x-dx/2, y, G, c_idx)/(dx*dx);
vec_jm(index-Nx) = -dt*D_func(x, y-dy/2, G, c_idx)/(dy*dy);

    
% Dirichlet boundary conditions for left and right boundary
% 3a) Set up rows for the left boundary
index = e_w;
vec(index) = 1;

% 4a) Set up rows for the right boundary
index = e_e;
vec(index) = 1;

% 1b) Set up rows for the north left boundary
index = e_nw;
vec(index) = 1;


% 2b) Set up rows for the north right boundary
index = e_ne;
vec(index) = 1;

% 3b) Set up rows for the south left boundary
index = e_sw;
vec(index) = 1;

% 4b) Set up rows for the south east boundary
index = e_se;
vec(index) = 1;

%%%%%%% INNER DOMAIN %%%%%%%
index = e;
[x, y] = x_y_from_idx(index, G);
vec(index) = 1 + dt*(D_func(x+dx/2, y, G, c_idx)/(dx*dx) + ...
    D_func(x-dx/2, y, G, c_idx)/(dx*dx) + ...
    D_func(x, y+dy/2, G, c_idx)/(dy*dy) + ...
    D_func(x, y-dy/2, G, c_idx)/(dy*dy));
vec_kp(index+1) = -dt*D_func(x+dx/2, y, G, c_idx)/(dx*dx);
vec_km(index-1) = -dt*D_func(x-dx/2, y, G, c_idx)/(dx*dx);
vec_jp(index+Nx) = -dt*D_func(x, y+dy/2, G, c_idx)/(dy*dy);
vec_jm(index-Nx) = -dt*D_func(x, y-dy/2, G, c_idx)/(dy*dy);

%%%%%%% SET UP THE MATRIX %%%%%%%
A = spdiags(vec, 0, N, N);
A = A + spdiags(vec_kp, 1, N, N);
A = A + spdiags(vec_km, -1, N, N);
A = A + spdiags(vec_jp, Nx, N, N);
A = A + spdiags(vec_jm, -Nx, N, N);

end

