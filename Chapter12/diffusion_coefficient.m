function D = diffusion_coefficient(x, y, G, c_idx)
%D = diffusion_coefficient(x, y, G, c_idx)

% Zero in O_m, D(c_idx) in O_i and O_e. Only nonzero for potassium in O_c.
D = (x < G.Li).*G.D(c_idx) + (x > G.Li+G.Lm).*G.D(c_idx) + ...
    (c_idx == 2).*(x >= G.Li).*(x <= G.Li+G.Lm).*(y > G.channel_start).*(y < G.channel_end).*G.D(c_idx);

end

