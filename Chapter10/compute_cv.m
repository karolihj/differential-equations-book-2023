function CV = compute_cv(v, v_th, t, x, y, x_start, x_end, y_start, y_end)
%CV = compute_cv(v, v_th, t, x, y, x_start, x_end, y_start, y_end)
% Compute the conduction velocity

% Find indices for measuring the cv
Mx = length(x);
[~, start_idx_x] = min(abs(x-x_start));
[~, end_idx_x] = min(abs(x-x_end));
[~, start_idx_y] = min(abs(y-y_start));
[~, end_idx_y] = min(abs(y-y_end));
start_idx = (start_idx_y-1)*Mx + start_idx_x;
end_idx = (end_idx_y-1)*Mx + end_idx_x;

t_start = 0;
t_end = inf;
cv_started = 0;
for n=1:length(t)

    % Find time point when v crosses threshold in x_start
    if ~cv_started && v(start_idx, n) >= v_th
        cv_started = 1;
        t_start = t(n);
    end

    % Find time point when v crosses threshold in x_end
    if v(end_idx, n) >= v_th
        t_end = t(n);
        break;
    end
end

% Compute conduction velocity
% Multiply by 1e-3 to convert from cm/ms to cm/s
CV = (sqrt((x(end_idx_x)-x(start_idx_x))^2+(y(end_idx_y)-y(start_idx_y))^2))/((t_end-t_start)*1e-3);

end