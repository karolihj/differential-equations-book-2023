function [x, y] = x_y_from_idx(idx, G)
%[x, y] = x_y_from_idx(idx, G)
% Get x and y value for the global index idx

x_idx = rem(idx-1, G.Nx);
y_idx = floor((idx-1)/G.Nx);

x = G.dx/2+G.dx*x_idx;
y = G.dy/2+G.dy*y_idx;
end

