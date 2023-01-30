function eps = permittivity(x, y, G)
%eps = permittivity(x, y, G)

% eps in O_i and O_e, eps_mem in O_m and O_c
eps = (x < G.Li).*G.eps + (x > G.Li+G.Lm)*G.eps + ...
    (x >= G.Li).*(x <= G.Li+G.Lm)*G.eps_mem;

end

