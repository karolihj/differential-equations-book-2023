function A = set_up_Cc_matrix(G, mesh, c_idx, phi)
%set_up_Cc_matrix Set up matrix for the finite difference equations for the
%concentrations (only the part that depends on phi)

phi = phi';

% Load parameters
N = G.N;
Nx = G.Nx;
elem = G.e;
kB = G.kB;
T = G.T;
z = G.z(c_idx);
dt = G.dt;
dx = G.dx;
dy = G.dy;
D_func = G.D_func;


% Load grid point names
e_s = mesh.e_s;
e_n = mesh.e_n;
e = mesh.e;

vec = zeros(N, 1);
vec_kp = zeros(N, 1);
vec_km = zeros(N, 1);
vec_jp = zeros(N, 1);
vec_jm = zeros(N, 1);


%%%%%%% BOUNDARY %%%%%%%   
% Neumann boundary condition for north and south boundaries
index = e_s;
[x, y] = x_y_from_idx(index, G);
vec(index) = (dt*D_func(x+dx/2, y, G, c_idx)*z*elem/(2*kB*T)).*(phi(index+1) - phi(index))/(dx*dx) + ...
    (dt*D_func(x-dx/2, y, G, c_idx)*z*elem/(2*kB*T)).*(phi(index-1) - phi(index))/(dx*dx) + ...
    (dt*D_func(x, y+dy/2, G, c_idx)*z*elem/(2*kB*T)).*(phi(index+Nx) - phi(index))/(dy*dy);
vec_kp(index+1) = (dt*D_func(x+dx/2, y, G, c_idx)*z*elem/(2*kB*T)).*(phi(index+1) - phi(index))/(dx*dx);
vec_km(index-1) = (dt*D_func(x-dx/2, y, G, c_idx)*z*elem/(2*kB*T)).*(phi(index-1) - phi(index))/(dx*dx);
vec_jp(index+Nx) = (dt*D_func(x, y+dy/2, G, c_idx)*z*elem/(2*kB*T)).*((phi(index+Nx) - phi(index))/(dy*dy));

index = e_n;
[x, y] = x_y_from_idx(index, G);
vec(index) = (dt*D_func(x+dx/2, y, G, c_idx)*z*elem/(2*kB*T)).*(phi(index+1) - phi(index))/(dx*dx) + ...
    (dt*D_func(x-dx/2, y, G, c_idx)*z*elem/(2*kB*T)).*(phi(index-1) - phi(index))/(dx*dx) + ...
    (dt*D_func(x, y-dy/2, G, c_idx)*z*elem/(2*kB*T)).*(phi(index-Nx) - phi(index))/(dy*dy);
vec_kp(index+1) = (dt*D_func(x+dx/2, y, G, c_idx)*z*elem/(2*kB*T)).*(phi(index+1) - phi(index))/(dx*dx);
vec_km(index-1) = (dt*D_func(x-dx/2, y, G, c_idx)*z*elem/(2*kB*T)).*(phi(index-1) - phi(index))/(dx*dx);
vec_jm(index-Nx) = (dt*D_func(x, y-dy/2, G, c_idx)*z*elem/(2*kB*T)).*((phi(index-Nx) - phi(index))/(dy*dy));


%%%%%%% INNER DOMAIN %%%%%%%
index = e;
[x, y] = x_y_from_idx(index, G);
vec(index) = (dt*D_func(x+dx/2, y, G, c_idx)*z*elem/(2*kB*T)).*(phi(index+1) - phi(index))/(dx*dx) + ...
    (dt*D_func(x-dx/2, y, G, c_idx)*z*elem/(2*kB*T)).*(phi(index-1) - phi(index))/(dx*dx) + ...
    (dt*D_func(x, y+dy/2, G, c_idx)*z*elem/(2*kB*T)).*(phi(index+Nx) - phi(index))/(dy*dy) + ...
    (dt*D_func(x, y-dy/2, G, c_idx)*z*elem/(2*kB*T)).*(phi(index-Nx) - phi(index))/(dy*dy);
vec_kp(index+1) = (dt*D_func(x+dx/2, y, G, c_idx)*z*elem/(2*kB*T)).*(phi(index+1) - phi(index))/(dx*dx);
vec_km(index-1) = (dt*D_func(x-dx/2, y, G, c_idx)*z*elem/(2*kB*T)).*(phi(index-1) - phi(index))/(dx*dx);
vec_jp(index+Nx) = (dt*D_func(x, y+dy/2, G, c_idx)*z*elem/(2*kB*T)).*((phi(index+Nx) - phi(index))/(dy*dy));
vec_jm(index-Nx) = (dt*D_func(x, y-dy/2, G, c_idx)*z*elem/(2*kB*T)).*((phi(index-Nx) - phi(index))/(dy*dy));

%%%%%%% SET UP THE MATRIX %%%%%%%
A = spdiags(vec, 0, N, N);
A = A + spdiags(vec_kp, 1, N, N);
A = A + spdiags(vec_km, -1, N, N);
A = A + spdiags(vec_jp, Nx, N, N);
A = A + spdiags(vec_jm, -Nx, N, N);

A = -A;

end

