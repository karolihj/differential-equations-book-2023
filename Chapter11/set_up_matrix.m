function A = set_up_matrix(G, mesh)
%set_up_matrix Set up matrix for the finite difference equations

% Load parameters
N = G.N;
Nx = G.Nx;
Ny = G.Ny;
dx = G.dx;
dy = G.dy;
dz= G.dz;
Cm = G.Cm;
dt = G.dt;
sigma_e = G.sigma_e;
sigma_i = G.sigma_i;

% Load grid point names
e_lsw = mesh.e_lsw;
e_lse = mesh.e_lse;
e_lnw = mesh.e_lnw;
e_lne = mesh.e_lne;
e_hsw = mesh.e_hsw;
e_hse = mesh.e_hse;
e_hnw = mesh.e_hnw;
e_hne = mesh.e_hne;

e_hw = mesh.e_hw;
e_he = mesh.e_he;
e_hs = mesh.e_hs;
e_hn = mesh.e_hn;
e_lw = mesh.e_lw;
e_le = mesh.e_le;
e_ls = mesh.e_ls;
e_ln = mesh.e_ln;
e_ne = mesh.e_ne;
e_sw = mesh.e_sw;
e_se = mesh.e_se;
e_nw = mesh.e_nw;

e_w = mesh.e_w;
e_e = mesh.e_e;
e_s = mesh.e_s;
e_n = mesh.e_n;
e_h = mesh.e_h;
e_l = mesh.e_l;

m_lsw = mesh.m_lsw;
m_lse = mesh.m_lse;
m_lnw = mesh.m_lnw;
m_lne = mesh.m_lne;
m_hsw = mesh.m_hsw;
m_hse = mesh.m_hse;
m_hnw = mesh.m_hnw;
m_hne = mesh.m_hne;

m_hw = mesh.m_hw;
m_he = mesh.m_he;
m_hs = mesh.m_hs;
m_hn = mesh.m_hn;
m_lw = mesh.m_lw;
m_le = mesh.m_le;
m_ls = mesh.m_ls;
m_ln = mesh.m_ln;
m_ne = mesh.m_ne;
m_sw = mesh.m_sw;
m_se = mesh.m_se;
m_nw = mesh.m_nw;

m_w = mesh.m_w;
m_e = mesh.m_e;
m_s = mesh.m_s;
m_n = mesh.m_n;
m_h = mesh.m_h;
m_l = mesh.m_l;

i_lsw = mesh.i_lsw;
i_lse = mesh.i_lse;
i_lnw = mesh.i_lnw;
i_lne = mesh.i_lne;
i_hsw = mesh.i_hsw;
i_hse = mesh.i_hse;
i_hnw = mesh.i_hnw;
i_hne = mesh.i_hne;

i_hw = mesh.i_hw;
i_he = mesh.i_he;
i_hs = mesh.i_hs;
i_hn = mesh.i_hn;
i_lw = mesh.i_lw;
i_le = mesh.i_le;
i_ls = mesh.i_ls;
i_ln = mesh.i_ln;
i_ne = mesh.i_ne;
i_sw = mesh.i_sw;
i_se = mesh.i_se;
i_nw = mesh.i_nw;

i_w = mesh.i_w;
i_e = mesh.i_e;
i_s = mesh.i_s;
i_n = mesh.i_n;
i_h = mesh.i_h;
i_l = mesh.i_l;

i = mesh.i;
e = mesh.e;


vec = zeros(2*N, 1);
vec_kp = zeros(2*N, 1);
vec_km = zeros(2*N, 1);
vec_jp = zeros(2*N, 1);
vec_jm = zeros(2*N, 1);
vec_qp = zeros(2*N, 1);
vec_qm = zeros(2*N, 1);

vec_v = zeros(2*N, 1);
vec_km_v = zeros(2*N, 1);
vec_kp_v = zeros(2*N, 1);
vec_jm_v = zeros(2*N, 1);
vec_jp_v = zeros(2*N, 1);
vec_qm_v = zeros(2*N, 1);
vec_qp_v = zeros(2*N, 1);

v_vec = zeros(2*N, 1);
v_vec_kp = zeros(2*N, 1);
v_vec_km = zeros(2*N, 1);
v_vec_jp = zeros(2*N, 1);
v_vec_jm = zeros(2*N, 1);
v_vec_qp = zeros(2*N, 1);
v_vec_qm = zeros(2*N, 1);


% Set up rows for the intracellular domain
index = [i, i_lsw, i_lse, i_lnw, i_lne, i_hsw, i_hse, i_hnw, i_hne, i_hw, i_he i_hs, i_hn, ...
    i_lw, i_le, i_ls, i_ln, i_ne, i_sw, i_se, i_nw, i_w, i_e, i_s, i_n, i_h, i_l];
vec(index) = -((sigma_i+sigma_i)/(dx*dx) + (sigma_i+sigma_i)/(dy*dy) + (sigma_i+sigma_i)/(dz*dz)); 
vec_kp(index+1) = sigma_i/(dx*dx);
vec_km(index-1) = sigma_i/(dx*dx);
vec_jp(index+Nx) = sigma_i/(dy*dy);
vec_jm(index-Nx) = sigma_i/(dy*dy);
vec_qp(index+Nx*Ny) = sigma_i/(dz*dz);
vec_qm(index-Nx*Ny) = sigma_i/(dz*dz);


% Set up extra term for the high inner boundary
index = [i_h, i_hw, i_he, i_hs, i_hn, i_hsw, i_hse, i_hnw, i_hne];
vec_qp_v(N+index+Nx*Ny) = sigma_i/(dz*dz);

% Set up extra term for the low inner boundary
index = [i_l, i_lw, i_le, i_ls, i_ln, i_lsw, i_lse, i_lnw, i_lne];
vec_qm_v(N+index-Nx*Ny) = sigma_i/(dz*dz);

% Set up extra term for the north inner boundary
index = [i_n, i_nw, i_ne, i_hn, i_ln, i_lnw, i_lne, i_hnw, i_hne];
vec_jp_v(N+index+Nx) = sigma_i/(dy*dy);

% Set up extra term for the south inner boundary
index = [i_s, i_sw, i_se, i_hs, i_ls, i_lsw, i_lse, i_hsw, i_hse];
vec_jm_v(N+index-Nx) = sigma_i/(dy*dy);

% Set up extra term for the right inner boundary
index = [i_e, i_ne, i_se, i_he, i_le, i_lne, i_lse, i_hse, i_hne];
vec_kp_v(N+index+1) = sigma_i/(dx*dx);

% Set up extra term for the left inner boundary
index = [i_w, i_nw, i_sw, i_hw, i_lw, i_lnw, i_lsw, i_hsw, i_hnw];
vec_km_v(N+index-1) = sigma_i/(dx*dx);


% Set up terms for the high membrane boundary
index = [m_h, m_hw, m_he, m_hs, m_hn, m_hsw, m_hse, m_hnw, m_hne];
vec(index) = vec(index) -((sigma_e+sigma_i)/dz);
vec_qp(index+Nx*Ny) = vec_qp(index+Nx*Ny) + sigma_e/dz;
vec_qm(index-Nx*Ny) = vec_qm(index-Nx*Ny) + sigma_i/dz;
vec_v(N+index) = vec_v(N+index) - sigma_i/dz;
vec_qm_v(N+index-Nx*Ny) = vec_qm_v(N+index-Nx*Ny) + sigma_i/dz;

% Set up terms for the low membrane boundary
index = [m_l, m_lw, m_le, m_ls, m_ln, m_lsw, m_lse, m_lnw, m_lne];
vec(index) = vec(index) - ((sigma_e+sigma_i)/dz);
vec_qp(index+Nx*Ny) = vec_qp(index+Nx*Ny) + sigma_i/dz;
vec_qm(index-Nx*Ny) = vec_qm(index-Nx*Ny) + sigma_e/dz;
vec_v(N+index) = vec_v(N+index) - sigma_i/dz;
vec_qp_v(N+index+Nx*Ny) = vec_qp_v(N+index+Nx*Ny) + sigma_i/dz;

% Set up terms for the north membrane boundary
index = [m_n, m_nw, m_ne, m_hn, m_ln, m_lnw, m_lne, m_hnw, m_hne];
vec(index) = vec(index) -((sigma_e+sigma_i)/dy);
vec_jp(index+Nx) = vec_jp(index+Nx) + sigma_e/dy;
vec_jm(index-Nx) = vec_jm(index-Nx) + sigma_i/dy;
vec_v(N+index) = vec_v(N+index) -sigma_i/dy;
vec_jm_v(N+index-Nx) = vec_jm_v(N+index-Nx) + sigma_i/dy;


% Set up terms for the south membrane boundary
index = [m_s, m_sw, m_se, m_hs, m_ls, m_lsw, m_lse, m_hsw, m_hse];
vec(index) = vec(index) -(sigma_e+sigma_i)/dy;
vec_jp(index+Nx) = vec_jp(index+Nx) + sigma_i/dy;
vec_jm(index-Nx) = vec_jm(index-Nx) + sigma_e/dy;
vec_v(N+index) = vec_v(N+index) - sigma_i/dy;
vec_jp_v(N+index+Nx) = vec_jp_v(N+index+Nx) + sigma_i/dy;


% Set up terms for the right membrane boundary
index = [m_e, m_ne, m_se, m_he, m_le, m_lne, m_lse, m_hse, m_hne];
vec(index) = vec(index) -((sigma_e+sigma_i)/dx);
vec_kp(index+1) = vec_kp(index+1) + sigma_e/dx;
vec_km(index-1) = vec_km(index-1) + sigma_i/dx;
vec_v(N+index) = vec_v(N+index) -sigma_i/dx;
vec_km_v(N+index-1) = vec_km_v(N+index-1) + sigma_i/dx;

% Set up terms for the left membrane boundary
index = [m_w, m_nw, m_sw, m_hw, m_lw, m_lnw, m_lsw, m_hsw, m_hnw];
vec(index) = vec(index) -((sigma_i+sigma_e)/dx);
vec_kp(index+1) = vec_kp(index+1) + sigma_i/dx;
vec_km(index-1) = vec_km(index-1) + sigma_e/dx;
vec_v(N+index) = vec_v(N+index) - sigma_i/dx;
vec_kp_v(N+index+1) = vec_kp_v(N+index+1) + sigma_i/dx;


% Extracellular boundary
index = [e_w, e_e, e_n, e_s, e_h, e_l, e_hw, e_he, e_hs, e_hn, e_lw, ...
    e_le, e_ls, e_ln, e_nw, e_ne, e_sw, e_se, e_lsw, e_lse, e_lnw, ...
    e_lne, e_hsw, e_hse, e_hnw, e_hne];
vec(index) = 1;
 
% Set up rows for the inner extracellular domain 
index = e;
vec(index) = -((sigma_e+sigma_e)/(dx*dx) + (sigma_e+sigma_e)/(dy*dy) ...
                    + (sigma_e+sigma_e)/(dz*dz));
vec_kp(index+1) = sigma_e/(dx*dx);
vec_km(index-1) = sigma_e/(dx*dx);
vec_jp(index+Nx) = sigma_e/(dy*dy);
vec_jm(index-Nx) = sigma_e/(dy*dy);
vec_qp(index+Nx*Ny) = sigma_e/(dz*dz);
vec_qm(index-Nx*Ny) = sigma_e/(dz*dz);

% Set up rows for the membrane potential
index = [m_w, m_e, m_n, m_s, m_h, m_l, m_hw, m_he, m_hs, m_hn, ...
    m_lw, m_le, m_ls, m_ln, m_nw, m_ne, m_sw, m_se, m_lsw, m_lse, ...
    m_lnw, m_lne, m_hsw, m_hse, m_hnw, m_hne];
vec(N+index) = Cm/dt;

% 1a) Set up factors for the low membrane
index = m_l;
v_vec(index) = -sigma_e/dz;
v_vec_qm(index-Nx*Ny) = sigma_e/dz;

% 2a) Set up factors for the high membrane
index = m_h;
v_vec(index) = -sigma_e/dz;
v_vec_qp(index+Nx*Ny) = sigma_e/dz;

% 3a) Set up factors for the south membrane
index = m_s;
v_vec(index) = -sigma_e/dy;
v_vec_jm(index-Nx) = sigma_e/dy;

% 4a) Set up factors for the north membrane
index = m_n;
v_vec(index) = -sigma_e/dy;
v_vec_jp(index+Nx) = sigma_e/dy;

% 5a) Set up factors for the left membrane
index = m_w;
v_vec(index) = -sigma_e/dx;
v_vec_km(index-1) = sigma_e/dx;

% 6a) Set up factors for the right membrane
index = m_e;
v_vec(index) = -sigma_e/dx;
v_vec_kp(index+1) = sigma_e/dx;

% 1b) Set up factors for the high left membrane
index = m_hw;
v_vec(index) = -0.5*(sigma_e/dz+sigma_e/dx);
v_vec_qp(index+Nx*Ny) = 0.5*sigma_e/dz;
v_vec_km(index-1) = 0.5*sigma_e/dx;

% 2b) Set up factors for the high right membrane
index = m_he;
v_vec(index) = -0.5*(sigma_e/dz+sigma_e/dx);
v_vec_qp(index+Nx*Ny) = 0.5*sigma_e/dz;
v_vec_kp(index+1) = 0.5*sigma_e/dx;

% 3b) Set up factors for the high south membrane
index = m_hs;
v_vec(index) = -0.5*(sigma_e/dz+sigma_e/dy);
v_vec_qp(index+Nx*Ny) = 0.5*sigma_e/dz;
v_vec_jm(index-Nx) = 0.5*sigma_e/dy;

% 4b) Set up factors for the high north membrane
index = m_hn;
v_vec(index) = -0.5*(sigma_e/dz+sigma_e/dy);
v_vec_qp(index+Nx*Ny) = 0.5*sigma_e/dz;
v_vec_jp(index+Nx) = 0.5*sigma_e/dy;

% 5b) Set up factors for the low left membrane
index = m_lw;
v_vec(index) = -0.5*(sigma_e/dz+sigma_e/dx);
v_vec_qm(index-Nx*Ny) = 0.5*sigma_e/dz;
v_vec_km(index-1) = 0.5*sigma_e/dx;

% 6b) Set up factors for the low right membrane
index = m_le;
v_vec(index) = -0.5*(sigma_e/dz+sigma_e/dx);
v_vec_qm(index-Nx*Ny) = 0.5*sigma_e/dz;
v_vec_kp(index+1) = 0.5*sigma_e/dx;

% 7b) Set up factors for the low south membrane
index = m_ls;
v_vec(index) = -0.5*(sigma_e/dz+sigma_e/dy);
v_vec_qm(index-Nx*Ny) = 0.5*sigma_e/dz;
v_vec_jm(index-Nx) = 0.5*sigma_e/dy;

% 8b) Set up factors for the low north membrane
index = m_ln;
v_vec(index) = -0.5*(sigma_e/dz+sigma_e/dy);
v_vec_qm(index-Nx*Ny) = 0.5*sigma_e/dz;
v_vec_jp(index+Nx) = 0.5*sigma_e/dy;

% 9b) Set up factors for the north left membrane
index = m_nw;
v_vec(index) = -0.5*(sigma_e/dy+sigma_e/dx);
v_vec_jp(index+Nx) = 0.5*sigma_e/dy;
v_vec_km(index-1) = 0.5*sigma_e/dx;

% 10b) Set up factors for the north right membrane
index = m_ne;
v_vec(index) = -0.5*(sigma_e/dy+sigma_e/dx);
v_vec_jp(index+Nx) = 0.5*sigma_e/dy;
v_vec_kp(index+1) = 0.5*sigma_e/dx;

% 11b) Set up factors for the south left membrane
index = m_sw;
v_vec(index) = -0.5*(sigma_e/dy+sigma_e/dx);
v_vec_jm(index-Nx) = 0.5*sigma_e/dy;
v_vec_km(index-1) = 0.5*sigma_e/dx;

% 12b) Set up factors for the south east membrane
index = m_se;
v_vec(index) = -0.5*(sigma_e/dy+sigma_e/dx);
v_vec_jm(index-Nx) = 0.5*sigma_e/dy;
v_vec_kp(index+1) = 0.5*sigma_e/dx;

% 1c) Set up factors for the lower, south, left membrane
index = m_lsw;
v_vec(index) = -1/3*(sigma_e/dz+sigma_e/dy+sigma_e/dx);
v_vec_qm(index-Nx*Ny) = 1/3*sigma_e/dz;
v_vec_jm(index-Nx) = 1/3*sigma_e/dy;
v_vec_km(index-1) = 1/3*sigma_e/dx;

% 2c) Set up factors for the lower, south, east membrane
index = m_lse;
v_vec(index) = -1/3*(sigma_e/dz+sigma_e/dy+sigma_e/dx);
v_vec_qm(index-Nx*Ny) = 1/3*sigma_e/dz;
v_vec_jm(index-Nx) = 1/3*sigma_e/dy;
v_vec_kp(index+1) = 1/3*sigma_e/dx;

% 3c) Set up factors for the lower, north, left membrane
index = m_lnw;
v_vec(index) = -1/3*(sigma_e/dz+sigma_e/dy+sigma_e/dx);
v_vec_qm(index-Nx*Ny) = 1/3*sigma_e/dz;
v_vec_jp(index+Nx) = 1/3*sigma_e/dy;
v_vec_km(index-1) = 1/3*sigma_e/dx;

% 4c) Set up factors for the lower, north, right membrane
index = m_lne;
v_vec(index) = -1/3*(sigma_e/dz+sigma_e/dy+sigma_e/dx);
v_vec_qm(index-Nx*Ny) = 1/3*sigma_e/dz;
v_vec_jp(index+Nx) = 1/3*sigma_e/dy;
v_vec_kp(index+1) = 1/3*sigma_e/dx;

% 5c) Set up factors for the higher, south, left membrane
index = m_hsw;
v_vec(index) = -1/3*(sigma_e/dz+sigma_e/dy+sigma_e/dx);
v_vec_qp(index+Nx*Ny) = 1/3*sigma_e/dz;
v_vec_jm(index-Nx) = 1/3*sigma_e/dy;
v_vec_km(index-1) = 1/3*sigma_e/dx;

% 6c) Set up factors for the higher, south, east membrane
index = m_hse;
v_vec(index) = -1/3*(sigma_e/dz+sigma_e/dy+sigma_e/dx);
v_vec_qp(index+Nx*Ny) = 1/3*sigma_e/dz;
v_vec_jm(index-Nx) = 1/3*sigma_e/dy;
v_vec_kp(index+1) = 1/3*sigma_e/dx;

% 7c) Set up factors for the higher, north, left membrane
index = m_hnw;
v_vec(index) = -1/3*(sigma_e/dz+sigma_e/dy+sigma_e/dx);
v_vec_qp(index+Nx*Ny) = 1/3*sigma_e/dz;
v_vec_jp(index+Nx) = 1/3*sigma_e/dy;
v_vec_km(index-1) = 1/3*sigma_e/dx;

% 8c) Set up factors for the higher, north, right membrane
index = m_hne;
v_vec(index) = -1/3*(sigma_e/dz+sigma_e/dy+sigma_e/dx);
v_vec_qp(index+Nx*Ny) = 1/3*sigma_e/dz;
v_vec_jp(index+Nx) = 1/3*sigma_e/dy;
v_vec_kp(index+1) = 1/3*sigma_e/dx;


% Set up the matrix
A = spdiags(vec, 0, 2*N, 2*N);
A = A + spdiags(vec_kp, 1, 2*N, 2*N);
A = A + spdiags(vec_km, -1, 2*N, 2*N);
A = A + spdiags(vec_jp, Nx, 2*N, 2*N);
A = A + spdiags(vec_jm, -Nx, 2*N, 2*N);
A = A + spdiags(vec_qp, Nx*Ny, 2*N, 2*N);
A = A + spdiags(vec_qm, -Nx*Ny, 2*N, 2*N);

A = A + spdiags(vec_v, N, 2*N, 2*N);
A = A + spdiags(vec_km_v, N-1, 2*N, 2*N);
A = A + spdiags(vec_kp_v, N+1, 2*N, 2*N);
A = A + spdiags(vec_jp_v, N+Nx, 2*N, 2*N);
A = A + spdiags(vec_jm_v, N-Nx, 2*N, 2*N);
A = A + spdiags(vec_qp_v, N+Nx*Ny, 2*N, 2*N);
A = A + spdiags(vec_qm_v, N-Nx*Ny, 2*N, 2*N);

A = A + spdiags(v_vec, -N, 2*N, 2*N);
A = A + spdiags(v_vec_km, -N-1, 2*N, 2*N);
A = A + spdiags(v_vec_kp, -N+1, 2*N, 2*N);
A = A + spdiags(v_vec_jm, -N-Nx, 2*N, 2*N);
A = A + spdiags(v_vec_jp, -N+Nx, 2*N, 2*N);
A = A + spdiags(v_vec_qm, -N-Nx*Ny, 2*N, 2*N);
A = A + spdiags(v_vec_qp, -N+Nx*Ny, 2*N, 2*N);



% Reshape matrix to get fewer unknowns
A_top = A(1:N, [1:N, N+mesh.v]);
A_bottom = A(N+1:2*N, [1:N, N+mesh.v]);
A_bottom = A_bottom(mesh.v,:);
A = [A_top; A_bottom];

end

