function CV = compute_cv(v, v_th, t, x, x_start, x_end)
%CV = compute_cv(v, v_th, t, x, x_start, x_end)
% Compute the conduction velocity

% Find indices for measuring the cv
[~, start_idx] = min(abs(x-x_start));
[~, end_idx] = min(abs(x-x_end));

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

CV = (x(end_idx)-x(start_idx))/((t_end-t_start)*1e-3);

end
