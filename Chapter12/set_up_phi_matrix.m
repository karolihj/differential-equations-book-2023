function A = set_up_phi_matrix(G, mesh)
%set_up_matrix Set up matrix for the finite difference equation for phi

% Load parameters
N = G.N;
Nx = G.Nx;
dx = G.dx;
dy = G.dy;
eps_func = G.eps_func;

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

% Dirichlet bc in east
index = [e_e, e_ne, e_se];
vec(index) = 1;

% Neumann on west, north and south boundaries
index = e_w;
[x, y] = x_y_from_idx(index, G);
vec(index) = -(eps_func(x+dx/2, y, G)/(dx*dx) + eps_func(x+dx/2, y, G)/(dx*dx) ...
+ eps_func(x, y+dy/2, G)/(dy*dy) + eps_func(x, y-dy/2, G)/(dy*dy)); 
vec_kp(index+1) = 2*eps_func(x+dx/2, y, G)/(dx*dx);
vec_jp(index+Nx) = eps_func(x, y+dy/2, G)/(dy*dy);
vec_jm(index-Nx) = eps_func(x, y-dy/2, G)/(dy*dy);

index = e_s;
[x, y] = x_y_from_idx(index, G);
vec(index) = -(eps_func(x+dx/2, y, G)/(dx*dx) + eps_func(x-dx/2, y, G)/(dx*dx) ...
    + eps_func(x, y+dy/2, G)/(dy*dy) + eps_func(x, y+dy/2, G)/(dy*dy));
vec_kp(index+1) = eps_func(x+dx/2, y, G)/(dx*dx);
vec_km(index-1) = eps_func(x-dx/2, y, G)/(dx*dx);
vec_jp(index+Nx) = 2*eps_func(x, y+dy/2, G)/(dy*dy);

index = e_n;
[x, y] = x_y_from_idx(index, G);
vec(index) = -(eps_func(x+dx/2, y, G)/(dx*dx) + eps_func(x-dx/2, y, G)/(dx*dx) ...
    + eps_func(x, y-dy/2, G)/(dy*dy) + eps_func(x, y-dy/2, G)/(dy*dy));
vec_kp(index+1) = eps_func(x+dx/2, y, G)/(dx*dx);
vec_km(index-1) = eps_func(x-dx/2, y, G)/(dx*dx);
vec_jm(index-Nx) = 2*eps_func(x, y-dy/2, G)/(dy*dy);

index = e_sw;
[x, y] = x_y_from_idx(index, G);
vec(index) = -(eps_func(x+dx/2, y, G)/(dx*dx) + eps_func(x+dx/2, y, G)/(dx*dx) ...
    + eps_func(x, y+dy/2, G)/(dy*dy) + eps_func(x, y+dy/2, G)/(dy*dy));
vec_kp(index+1) = 2*eps_func(x+dx/2, y, G)/(dx*dx);
vec_jp(index+Nx) = 2*eps_func(x, y+dy/2, G)/(dy*dy);

index = e_nw;
[x, y] = x_y_from_idx(index, G);
vec(index) = -(eps_func(x+dx/2, y, G)/(dx*dx) + eps_func(x+dx/2, y, G)/(dx*dx) ...
    + eps_func(x, y-dy/2, G)/(dy*dy) + eps_func(x, y-dy/2, G)/(dy*dy));
vec_kp(index+1) = 2*eps_func(x+dx/2, y, G)/(dx*dx);
vec_jm(index-Nx) = 2*eps_func(x, y-dy/2, G)/(dy*dy);


%%%%%%% INNER DOMAIN %%%%%%%
index = e;
[x, y] = x_y_from_idx(index, G);
vec(index) = -(eps_func(x+dx/2, y, G)/(dx*dx) + eps_func(x-dx/2, y, G)/(dx*dx) ...
    + eps_func(x, y+dy/2, G)/(dy*dy) + eps_func(x, y-dy/2, G)/(dy*dy));
vec_kp(index+1) = eps_func(x+dx/2, y, G)/(dx*dx);
vec_km(index-1) = eps_func(x-dx/2, y, G)/(dx*dx);
vec_jp(index+Nx) = eps_func(x, y+dy/2, G)/(dy*dy);
vec_jm(index-Nx) = eps_func(x, y-dy/2, G)/(dy*dy);

%%%%%%% SET UP THE MATRIX %%%%%%%
A = spdiags(vec, 0, N, N);
A = A + spdiags(vec_kp, 1, N, N);
A = A + spdiags(vec_km, -1, N, N);
A = A + spdiags(vec_jp, Nx, N, N);
A = A + spdiags(vec_jm, -Nx, N, N);

end

