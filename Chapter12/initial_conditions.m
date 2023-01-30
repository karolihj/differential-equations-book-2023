function c0 = initial_conditions(x, y, G, c_idx)
%c0 = initial_conditions(x, y, G, c_idx)

c0 = (x < G.Li)*G.c0i(c_idx) + (x > G.Li+G.Lm)*G.c0e(c_idx) + ...
    (c_idx==2)*(x >= G.Li).*(x <= G.Li+G.Lm).*...
    (y > G.channel_start-G.dy/2).*(y < G.channel_end+G.dy/2).*...
    (G.c0i(c_idx) + (x-(G.Li+G.dx/2))/(G.Lm-G.dx)*(G.c0e(c_idx)-G.c0i(c_idx))); 

end

