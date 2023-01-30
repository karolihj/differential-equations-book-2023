% Compute the numerical solution of the FitzHugh-Nagumo model with an added
% diffusion term. First using a fine time resolution, and then using two
% types of operator splitting
clear all

% Define model constants
a = -0.12;
c1 = 0.175;
c2 = 0.03;
b = 0.011;
d = 0.55;
delta = 5e-5;

% Define discretization parameters
L = 1;                % Length of spatial domain
dx = 0.01;            % Distance between spatial grid points
M = round(L/dx) + 1;  % Number of spatial grid points
x = (0:dx:L);         % x values


%% EXPLICIT SOLUTION WITH FINE TIME STEP, NO OPERATOR SPLITTING
dt = 1e-4;           % Time step
T = 100;             % Simulation time
N = round(T/dt);     % Number of time steps 
rho = delta*dt/dx^2;

% Set up arrays for saving the solutions
v_fine = zeros(M, N+1);
w_fine = zeros(M, N+1);

% Define initial conditions
v_fine(:, 1) = 0;
w_fine(:, 1) = 0;

% Define non-zero v at the left part of the domain
v_fine(1:5,1) = 0.26;

% Define diffusion matrix
D = spdiags((1-2*rho)*ones(M, 1), 0, M, M) + ...
    spdiags([0; 2*rho; rho*ones(M-2, 1)], 1, M, M) + ...
    spdiags([rho*ones(M-2, 1); 2*rho; 0], -1, M, M);


% Compute the numerical solution
for n=1:N
    v_fine(:,n+1) = D*v_fine(:, n) ...
        + dt*(c1*v_fine(:,n).*(v_fine(:,n)-a).*(1-v_fine(:,n)) - c2*w_fine(:,n));
    w_fine(:,n+1) = w_fine(:,n) + dt*(b*(v_fine(:,n)-d*w_fine(:,n)));

    % Print progress
    if rem(n-1,50000) == 0
        fprintf('%.0f%% done\n', 100*n/N/6)
    end
end

%% OPERATOR SPLITTING
dt_values = [5; 2; 1; 0.5; 0.2];
errors_first = zeros(length(dt_values), 1);
errors_second = zeros(length(dt_values), 1);


for i=1:length(dt_values)
    dt_op = dt_values(i);
    N = round(T/dt_op); 

    % Number of inner time steps in each operator splitting step
    N_per_op = round(dt_op/dt); 

    % Set up arrays for saving the solutions
    v_first = zeros(M, N+1);
    w_first = zeros(M, N+1);
    v_second = zeros(M, N+1);
    w_second = zeros(M, N+1);

    % Define initial conditions
    v_first(:,1) = 0;
    w_first(:,1) = 0;
    v_second(:,1) = 0;
    w_second(:,1) = 0;

    % Define non-zero v at the left part of the domain
    v_first(1:5,1) = 0.26;
    v_second(1:5,1) = 0.26;

    % Compute the numerical solution
    for n=1:N

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%    First order splitting    %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Step 1: Diffusion part
        for k=1:N_per_op
            v_first(:,n+1) = D*v_first(:,n);
            w_first(:,n+1) = w_first(:,n);
            v_first(:,n) = v_first(:,n+1); % Update for next iteration (k)
        end

        % Step 2: Reaction part
        for k=1:N_per_op
            v_first(:,n+1) = v_first(:,n+1) + ...
                dt*(c1*v_first(:,n+1).*(v_first(:,n+1)-a).*(1-v_first(:,n+1)) - c2*w_first(:,n+1));
            w_first(:,n+1) = w_first(:,n+1) + dt*(b*(v_first(:,n+1)-d*w_first(:,n+1)));
        end
    


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%    Second order splitting    %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Step 1: Diffusion step
        if n==1 % half step
            for k=1:N_per_op/2
                v_second(:,n+1) = D*v_second(:,n);
                w_second(:,n+1) = w_second(:,n);
                v_second(:,n) = v_second(:,n+1); % Update for next iteration (k)
            end
        else
            for k=1:N_per_op
                v_second(:,n+1) = D*v_second(:,n);
                w_second(:,n+1) = w_second(:,n);
                v_second(:,n) = v_second(:,n+1); % Update for next iteration (k)
            end
        end

        % Step 2: Reaction part
        for k=1:N_per_op
            v_second(:,n+1) = v_second(:,n+1) + ...
                dt*(c1*v_second(:,n+1).*(v_second(:,n+1)-a).*(1-v_second(:,n+1)) - c2*w_second(:,n+1));
            w_second(:,n+1) = w_second(:,n+1) + dt*(b*(v_second(:,n+1)-d*w_second(:,n+1)));
        end

        % Step 3: Half diffusion step for last time step
        if n==N
            for k=1:N_per_op/2
                v_second(:,n+1) = D*v_second(:,n+1);
            end
        end

        % Print progress
        if rem(n-1,round(5/dt_op)) == 0
            fprintf('%.0f%% done\n', i*100/6 + 100*n/N/6)
        end
    
    end

    % Compute error 
    errors_first(i) = max(abs(v_first(:,end)-v_fine(:,end)));
    errors_second(i) = max(abs(v_second(:,end)-v_fine(:,end)));

end


% Print the error
fprintf('\ndt      E1        E1/dt     E2        E2/dt^2\n')
for i=1:length(dt_values)
    fprintf('%-7g %-9.3g %-9.2g %-9.3g %-9.2g\n', dt_values(i), ...
        errors_first(i), errors_first(i)/dt_values(i), ...
        errors_second(i), errors_second(i)/dt_values(i)^2);
end


