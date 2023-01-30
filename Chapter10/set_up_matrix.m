function A = set_up_matrix(Mx, My, dx, dy, sigma_x, sigma_y)
%A = set_up_matrix(Mx, My, dx, dy, sigma_x, sigma_y)

% Define parameters
M = Mx*My;
rho_x = sigma_x/dx^2;
rho_y = sigma_y/dy^2;


% Define the matrix
A = zeros(M, M);
for j=1:My
    for k=1:Mx

        i = Mx*(j-1) + k; % Define global index

        if j==1 && k==1
            A(i, i) = -2*(rho_x + rho_y);
            A(i, i+1) = 2*rho_x;
            A(i, i+Mx) = 2*rho_y;
        elseif j==1 && k==Mx
            A(i, i-1) = 2*rho_x;
            A(i, i) = -2*(rho_x + rho_y);
            A(i, i+Mx) = 2*rho_y;
        elseif j==My && k==1
            A(i, i-Mx) = 2*rho_y;
            A(i, i) = -2*(rho_x + rho_y);
            A(i, i+1) = 2*rho_x;
        elseif j==My && k==Mx
            A(i, i-Mx) = 2*rho_y;
            A(i, i-1) = 2*rho_x;
            A(i, i) = -2*(rho_x + rho_y);
        elseif k==1
            A(i, i-Mx) = rho_y;
            A(i, i) = -2*(rho_x + rho_y);
            A(i, i+1) = 2*rho_x;
            A(i, i+Mx) = rho_y;
        elseif k==Mx
            A(i, i-Mx) = rho_y;
            A(i, i-1) = 2*rho_x;
            A(i, i) = -2*(rho_x + rho_y);
            A(i, i+Mx) = rho_y;
        elseif j==1
            A(i, i-1) = rho_x;
            A(i, i) = -2*(rho_x + rho_y);
            A(i, i+1) = rho_x;
            A(i, i+Mx) = 2*rho_y;
        elseif j==My
            A(i, i-Mx) = 2*rho_y;
            A(i, i-1) = rho_x;
            A(i, i) = -2*(rho_x + rho_y);
            A(i, i+1) = rho_x;
        else
            A(i, i-Mx) = rho_y;
            A(i, i-1) = rho_x;
            A(i, i) = -2*(rho_x + rho_y);
            A(i, i+1) = rho_x;
            A(i, i+Mx) = rho_y;
        end
    end
end

A = sparse(A);

end