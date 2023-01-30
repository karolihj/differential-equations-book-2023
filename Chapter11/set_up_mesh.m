function mesh = set_up_mesh(G)
% Set up vectors containing the indices of the different types of mesh nodes

% Read geometry
dx = G.dx;
dy = G.dy;
dz = G.dz;
Nx = G.Nx;
Ny = G.Ny;
Nz = G.Nz;
N = G.N;
cell_start_x = G.cell_start_x;
cell_start_y = G.cell_start_y;
cell_start_z = G.cell_start_z;
cell_length_x = G.cell_length_x;
cell_length_y = G.cell_length_y;
cell_length_z = G.cell_length_z;

% Define cell indices
cell_start_ix = round(cell_start_x/dx)+1;
cell_start_iy = round(cell_start_y/dy)+1;
cell_start_iz = round(cell_start_z/dz)+1;
cell_end_ix = cell_start_ix + round(cell_length_x/dx);
cell_end_iy = cell_start_iy + round(cell_length_y/dy);
cell_end_iz = cell_start_iz + round(cell_length_z/dz);

% Set up indices for outer boundary (corners)
mesh.e_lsw = 1;
mesh.e_lse = Nx;
mesh.e_lnw = (Ny-1)*Nx + 1;
mesh.e_lne = Nx*Ny;
mesh.e_hsw = (Nz-1)*Ny*Nx+1;
mesh.e_hse = (Nz-1)*Ny*Nx+Nx;
mesh.e_hnw = ((Nz-1)*Ny+Ny-1)*Nx + 1;
mesh.e_hne = Nx*Ny*Nz;

% Set up indices for outer boundary (lines)
mesh.e_lw = ((2:Ny-1)-1)*Nx + 1;
mesh.e_le = ((2:Ny-1)-1)*Nx + Nx;
mesh.e_ls = 2:Nx-1;
mesh.e_ln = (Ny-1)*Nx + (2:Nx-1);
mesh.e_hw = ((Nz-1)*Ny+(2:Ny-1)-1)*Nx + 1;
mesh.e_he = ((Nz-1)*Ny+(2:Ny-1)-1)*Nx + Nx;
mesh.e_hs = (Nz-1)*Ny*Nx+(2:Nx-1);
mesh.e_hn = ((Nz-1)*Ny+Ny-1)*Nx + (2:Nx-1);
mesh.e_sw = (1:Nz-2)*Ny*Nx+1;
mesh.e_se = (1:Nz-2)*Ny*Nx+Nx;
mesh.e_nw = ((1:Nz-2)*Ny+Ny-1)*Nx+1;
mesh.e_ne = ((1:Nz-2)*Ny+Ny-1)*Nx+Nx;

[x, y] = meshgrid(2:Nx-1, 2:Ny-1);
mesh.e_l = sort(reshape(sub2ind([Nx,Ny,Nz], x, y, ones((Ny-2),(Nx-2))), (Nx-2)*(Ny-2), 1))';
mesh.e_h = sort(reshape(sub2ind([Nx,Ny,Nz], x, y, Nz*ones((Ny-2),(Nx-2))), (Nx-2)*(Ny-2), 1))';
[x, z] = meshgrid(2:Nx-1, 2:Nz-1);
mesh.e_s = sort(reshape(sub2ind([Nx,Ny,Nz], x, ones((Nz-2), (Nx-2)), z), (Nx-2)*(Nz-2), 1))';
mesh.e_n = sort(reshape(sub2ind([Nx,Ny,Nz], x, Ny*ones((Nz-2), (Nx-2)), z), (Nx-2)*(Nz-2), 1))';
[y, z] = meshgrid(2:Ny-1, 2:Nz-1);
mesh.e_w = sort(reshape(sub2ind([Nx,Ny,Nz], ones((Nz-2), (Ny-2)), y, z), (Ny-2)*(Nz-2), 1))';
mesh.e_e =  sort(reshape(sub2ind([Nx,Ny,Nz], Nx*ones((Nz-2), (Ny-2)), y, z), (Ny-2)*(Nz-2), 1))';

% Set up indices for the membrane (corners)
mesh.m_lsw = ((cell_start_iz-1)*Ny+cell_start_iy-1)*Nx + cell_start_ix;
mesh.m_lse = ((cell_start_iz-1)*Ny+cell_start_iy-1)*Nx + cell_end_ix;
mesh.m_lnw = ((cell_start_iz-1)*Ny+cell_end_iy-1)*Nx + cell_start_ix;
mesh.m_lne = ((cell_start_iz-1)*Ny+cell_end_iy-1)*Nx + cell_end_ix;
mesh.m_hsw = ((cell_end_iz-1)*Ny+cell_start_iy-1)*Nx + cell_start_ix;
mesh.m_hse = ((cell_end_iz-1)*Ny+cell_start_iy-1)*Nx + cell_end_ix;
mesh.m_hnw = ((cell_end_iz-1)*Ny+cell_end_iy-1)*Nx + cell_start_ix;
mesh.m_hne = ((cell_end_iz-1)*Ny+cell_end_iy-1)*Nx + cell_end_ix;

% Set up indices for the membrane (lines)
mesh.m_ls = ((cell_start_iz-1)*Ny+cell_start_iy-1)*Nx + (cell_start_ix+1:cell_end_ix-1);
mesh.m_ln = ((cell_start_iz-1)*Ny+cell_end_iy-1)*Nx + (cell_start_ix+1:cell_end_ix-1);
mesh.m_lw = ((cell_start_iz-1)*Ny+(cell_start_iy+1:cell_end_iy-1)-1)*Nx + cell_start_ix;
mesh.m_le = ((cell_start_iz-1)*Ny+(cell_start_iy+1:cell_end_iy-1)-1)*Nx + cell_end_ix;
mesh.m_hs = ((cell_end_iz-1)*Ny+cell_start_iy-1)*Nx + (cell_start_ix+1:cell_end_ix-1);
mesh.m_hn = ((cell_end_iz-1)*Ny+cell_end_iy-1)*Nx + (cell_start_ix+1:cell_end_ix-1);
mesh.m_hw = ((cell_end_iz-1)*Ny+(cell_start_iy+1:cell_end_iy-1)-1)*Nx + cell_start_ix;
mesh.m_he = ((cell_end_iz-1)*Ny+(cell_start_iy+1:cell_end_iy-1)-1)*Nx + cell_end_ix;
mesh.m_sw = ((((cell_start_iz+1):(cell_end_iz-1))-1)*Ny+cell_start_iy-1)*Nx + cell_start_ix;
mesh.m_se = ((((cell_start_iz+1):(cell_end_iz-1))-1)*Ny+cell_start_iy-1)*Nx + cell_end_ix;
mesh.m_nw = ((((cell_start_iz+1):(cell_end_iz-1))-1)*Ny+cell_end_iy-1)*Nx + cell_start_ix;
mesh.m_ne = ((((cell_start_iz+1):(cell_end_iz-1))-1)*Ny+cell_end_iy-1)*Nx + cell_end_ix;

% Set up indices for the membrane (sides)
[x, y] = meshgrid(cell_start_ix+1:cell_end_ix-1, cell_start_iy+1:cell_end_iy-1);
mesh.m_l = sort(reshape(sub2ind([Nx,Ny,Nz], x, y, cell_start_iz*ones(size(x))), size(x,1)*size(x,2), 1))';
mesh.m_h = sort(reshape(sub2ind([Nx,Ny,Nz], x, y, cell_end_iz*ones(size(x))), size(x,1)*size(x,2), 1))';
[x, z] = meshgrid(cell_start_ix+1:cell_end_ix-1, cell_start_iz+1:cell_end_iz-1);
mesh.m_s = sort(reshape(sub2ind([Nx,Ny,Nz], x, cell_start_iy*ones(size(x)), z), size(x,1)*size(x,2), 1))';
mesh.m_n = sort(reshape(sub2ind([Nx,Ny,Nz], x, cell_end_iy*ones(size(x)), z), size(x,1)*size(x,2), 1))';
[y, z] = meshgrid(cell_start_iy+1:cell_end_iy-1, cell_start_iz+1:cell_end_iz-1);
mesh.m_w = sort(reshape(sub2ind([Nx,Ny,Nz], cell_start_ix*ones(size(y)), y, z), size(y,1)*size(y,2), 1))';
mesh.m_e =  sort(reshape(sub2ind([Nx,Ny,Nz], cell_end_ix*ones(size(y)), y, z), size(y,1)*size(y,2), 1))';


% Set up indices for the inner boundary of the cell (corners)
mesh.i_lsw = ((cell_start_iz)*Ny+cell_start_iy)*Nx + cell_start_ix+1;
mesh.i_lse = ((cell_start_iz)*Ny+cell_start_iy)*Nx + cell_end_ix-1;
mesh.i_lnw = ((cell_start_iz)*Ny+cell_end_iy-2)*Nx + cell_start_ix+1;
mesh.i_lne = ((cell_start_iz)*Ny+cell_end_iy-2)*Nx + cell_end_ix-1;
mesh.i_hsw = ((cell_end_iz-2)*Ny+cell_start_iy)*Nx + cell_start_ix+1;
mesh.i_hse = ((cell_end_iz-2)*Ny+cell_start_iy)*Nx + cell_end_ix-1;
mesh.i_hnw = ((cell_end_iz-2)*Ny+cell_end_iy-2)*Nx + cell_start_ix+1;
mesh.i_hne = ((cell_end_iz-2)*Ny+cell_end_iy-2)*Nx + cell_end_ix-1;

% Set up indices for the inner boundary of the cell (lines)
mesh.i_ls = ((cell_start_iz)*Ny+cell_start_iy)*Nx + (cell_start_ix+2:cell_end_ix-2);
mesh.i_ln = ((cell_start_iz)*Ny+cell_end_iy-2)*Nx + (cell_start_ix+2:cell_end_ix-2);
mesh.i_lw = ((cell_start_iz)*Ny+(cell_start_iy+2:cell_end_iy-2)-1)*Nx + cell_start_ix+1;
mesh.i_le = ((cell_start_iz)*Ny+(cell_start_iy+2:cell_end_iy-2)-1)*Nx + cell_end_ix-1;
mesh.i_hs = ((cell_end_iz-2)*Ny+cell_start_iy)*Nx + (cell_start_ix+2:cell_end_ix-2);
mesh.i_hn = ((cell_end_iz-2)*Ny+cell_end_iy-2)*Nx + (cell_start_ix+2:cell_end_ix-2);
mesh.i_hw = ((cell_end_iz-2)*Ny+(cell_start_iy+2:cell_end_iy-2)-1)*Nx + cell_start_ix+1;
mesh.i_he = ((cell_end_iz-2)*Ny+(cell_start_iy+2:cell_end_iy-2)-1)*Nx + cell_end_ix-1;
mesh.i_sw = ((((cell_start_iz+2):(cell_end_iz-2))-1)*Ny+cell_start_iy)*Nx + cell_start_ix+1;
mesh.i_se = ((((cell_start_iz+2):(cell_end_iz-2))-1)*Ny+cell_start_iy)*Nx + cell_end_ix-1;
mesh.i_nw = ((((cell_start_iz+2):(cell_end_iz-2))-1)*Ny+cell_end_iy-2)*Nx + cell_start_ix+1;
mesh.i_ne = ((((cell_start_iz+2):(cell_end_iz-2))-1)*Ny+cell_end_iy-2)*Nx + cell_end_ix-1;

% Set up indices for the inner boundary of the cell (sides)
[x, y] = meshgrid(cell_start_ix+2:cell_end_ix-2, cell_start_iy+2:cell_end_iy-2);
mesh.i_l = sort(reshape(sub2ind([Nx,Ny,Nz], x, y, (cell_start_iz+1)*ones(size(x))), size(x,1)*size(x,2), 1))';
mesh.i_h = sort(reshape(sub2ind([Nx,Ny,Nz], x, y, (cell_end_iz-1)*ones(size(x))), size(x,1)*size(x,2), 1))';
[x, z] = meshgrid(cell_start_ix+2:cell_end_ix-2, cell_start_iz+2:cell_end_iz-2);
mesh.i_s = sort(reshape(sub2ind([Nx,Ny,Nz], x, (cell_start_iy+1)*ones(size(x)), z), size(x,1)*size(x,2), 1))';
mesh.i_n = sort(reshape(sub2ind([Nx,Ny,Nz], x, (cell_end_iy-1)*ones(size(x)), z), size(x,1)*size(x,2), 1))';
[y, z] = meshgrid(cell_start_iy+2:cell_end_iy-2, cell_start_iz+2:cell_end_iz-2);
mesh.i_w = sort(reshape(sub2ind([Nx,Ny,Nz], (cell_start_ix+1)*ones(size(y)), y, z), size(y,1)*size(y,2), 1))';
mesh.i_e =  sort(reshape(sub2ind([Nx,Ny,Nz], (cell_end_ix-1)*ones(size(y)), y, z), size(y,1)*size(y,2), 1))';


% Set up indices for the remaining intracellular cells
[x, y, z] = meshgrid(cell_start_ix:cell_end_ix, cell_start_iy:cell_end_iy, ...
    cell_start_iz:cell_end_iz);
mesh.i_all = sort(reshape(sub2ind([Nx,Ny,Nz], x, y, z), size(x,1)*size(x,2)*size(x,3), 1));
indices = zeros(1, N);
indices(mesh.i_all) = 1;
indices([mesh.m_w, mesh.m_e, mesh.m_n, mesh.m_s, mesh.m_h, mesh.m_l, ...
    mesh.m_hw, mesh.m_he, mesh.m_hs, mesh.m_hn, mesh.m_lw, mesh.m_le, ...
    mesh.m_ls, mesh.m_ln, mesh.m_nw, mesh.m_ne, mesh.m_sw, mesh.m_se, ...
    mesh.m_lsw, mesh.m_lse, mesh.m_lnw, mesh.m_lne, ...
    mesh.m_hsw, mesh.m_hse, mesh.m_hnw, mesh.m_hne, ...
    mesh.i_w, mesh.i_e, mesh.i_n, mesh.i_s, mesh.i_h, mesh.i_l, ...
    mesh.i_hw, mesh.i_he, mesh.i_hs, mesh.i_hn, mesh.i_lw, mesh.i_le, ...
    mesh.i_ls, mesh.i_ln, mesh.i_nw, mesh.i_ne, mesh.i_sw, mesh.i_se, ...
    mesh.i_lsw, mesh.i_lse, mesh.i_lnw, mesh.i_lne, ...
    mesh.i_hsw, mesh.i_hse, mesh.i_hnw, mesh.i_hne]) = 0;
mesh.i = round(find(indices));


% Set up indices for the remaining extracellular cells
indices = ones(1, N);
indices([mesh.i_all', ...
    mesh.e_w, mesh.e_e, mesh.e_n, mesh.e_s, mesh.e_h, mesh.e_l, ...
    mesh.e_hw, mesh.e_he, mesh.e_hs, mesh.e_hn, mesh.e_lw, mesh.e_le, ...
    mesh.e_ls, mesh.e_ln, mesh.e_nw, mesh.e_ne, mesh.e_sw, mesh.e_se, ...
    mesh.e_lsw, mesh.e_lse, mesh.e_lnw, mesh.e_lne, ...
    mesh.e_hsw, mesh.e_hse, mesh.e_hnw, mesh.e_hne]) = 0;
mesh.e = round(find(indices));

% Set up indices for the membrane potential
mesh.v = sort([mesh.m_w, mesh.m_e, mesh.m_n, mesh.m_s, mesh.m_h, mesh.m_l, ...
    mesh.m_hw, mesh.m_he, mesh.m_hs, mesh.m_hn, mesh.m_lw, mesh.m_le, ...
    mesh.m_ls, mesh.m_ln, mesh.m_nw, mesh.m_ne, mesh.m_sw, mesh.m_se, ...
    mesh.m_lsw, mesh.m_lse, mesh.m_lnw, mesh.m_lne, ...
    mesh.m_hsw, mesh.m_hse, mesh.m_hnw, mesh.m_hne]);

% Set up indices for all extracellular points
mesh.e_all = sort([mesh.e, mesh.e_w, mesh.e_e, mesh.e_n, mesh.e_s, mesh.e_h, mesh.e_l, ...
    mesh.e_hw, mesh.e_he, mesh.e_hs, mesh.e_hn, mesh.e_lw, mesh.e_le, ...
    mesh.e_ls, mesh.e_ln, mesh.e_nw, mesh.e_ne, mesh.e_sw, mesh.e_se, ...
    mesh.e_lsw, mesh.e_lse, mesh.e_lnw, mesh.e_lne, ...
    mesh.e_hsw, mesh.e_hse, mesh.e_hnw, mesh.e_hne, ...
    mesh.m_w, mesh.m_e, mesh.m_n, mesh.m_s, mesh.m_h, mesh.m_l, ...
    mesh.m_hw, mesh.m_he, mesh.m_hs, mesh.m_hn, mesh.m_lw, mesh.m_le, ...
    mesh.m_ls, mesh.m_ln, mesh.m_nw, mesh.m_ne, mesh.m_sw, mesh.m_se, ...
    mesh.m_lsw, mesh.m_lse, mesh.m_lnw, mesh.m_lne, ...
    mesh.m_hsw, mesh.m_hse, mesh.m_hnw, mesh.m_hne]);


% Set up indices for all nodes except the outer boundary
indices = ones(1,N);
indices([mesh.e_w, mesh.e_e, mesh.e_n, mesh.e_s, mesh.e_h, mesh.e_l, ...
    mesh.e_hw, mesh.e_he, mesh.e_hs, mesh.e_hn, mesh.e_lw, mesh.e_le, ...
    mesh.e_ls, mesh.e_ln, mesh.e_nw, mesh.e_ne, mesh.e_sw, mesh.e_se, ...
    mesh.e_lsw, mesh.e_lse, mesh.e_lnw, mesh.e_lne, ...
    mesh.e_hsw, mesh.e_hse, mesh.e_hnw, mesh.e_hne]) = 0;
mesh.inner = round(find(indices));

% Set up indices to stimulate
num_stim = G.stim_length/G.dx;
[x, y, z] = meshgrid(cell_start_ix:(cell_start_ix+num_stim-1), cell_start_iy:cell_end_iy, ...
    cell_start_iz:cell_end_iz);
to_stim = sort(reshape(sub2ind([Nx,Ny,Nz], x, y, z), size(x,1)*size(x,2)*size(x,3), 1));
indices = zeros(1, N);
indices(to_stim) = indices(to_stim) + 1;
indices = indices(mesh.v);

mesh.to_stim = round(find(indices));



end