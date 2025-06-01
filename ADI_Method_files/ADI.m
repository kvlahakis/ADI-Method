function U = ADI(U_init, g, h, k, t_N)

%   Computes the finite‐difference approximation
%   of
%       u_t = u_{xx} + u_{yy},   (x,y) ∈ [0,1]×[0,1],   0 ≤ t ≤ t_N,
%   subject to an initial interior field U_init and Dirichlet boundary data 
%   g using The Peaceman–Rachford ADI method. A uniform discretization of 
%   the spatial domain is used giving {x_i, y_j} for 0 ≤ i,j ≤ m

%   Inputs
%       U_init  —  (m-1) × (m-1) matrix of initial interior values. U_init
%                  must be a matrix with entries (U_init)_ij = u(x_i,y_j,0)
%                  where 1 ≤ i,j ≤ m-1 and where m = floor(1/h).
%
%       g       —  function hand specifying the Dirichlet boundary conditions

%
%       h       —  Spatial step size in both x and y.
%
%       k       —  Temporal step size.
%
%       t_N     —  Final time.
%
%   Outputs
%       U       —  (m+1) x (m+1) matrix of numerical solution 
%                  on the entire grid at the final time. U is a matrix
%                  with entries U_ij ≈ u(x_i-1,y_j-1,t_N)

    U = U_init;
    m = floor(1/h);
    x = (0:h:1).';
    y = x;
    
    % operates on entire grid
    I = eye(m+1);
    e = ones(1, m);
    A =  1/(h^2) * (-2 * I + diag(e, -1) + diag(e, 1));
    M_plus = (I + k/2 * A);
    M_minus = (I - k/2 * A);

    % operates only on interior points
    I_tilde = I(2:end-1, 2:end-1);
    A_tilde = A(2:end-1, 2:end-1); 
    M_pl_tilde = (I_tilde + k/2 * A_tilde);
    M_min_tilde = (I_tilde - k/2 * A_tilde);

    t = 0:k:t_N;
    n = length(t);
    for i = 1:n-1
        t_n = t(i);
        t_np1 = t_n + k;
        
        % boundary condition for u star
        G_star = zeros(m-1, m-1);
        g_star_x_0 = 1/2 * M_plus * g(0, y, t_n) + 1/2 * M_minus * g(0, y, t_np1);
        g_star_x_1 = 1/2 * M_plus *g(1, y, t_n) + 1/2 * M_minus *g(1, y, t_np1);
        G_star(1,:) = g_star_x_0(2:end-1);
        G_star(end,:) = g_star_x_1(2:end-1);
    
        % boundary conditions (only interior points needed)
        G_n_y = zeros(m-1, m-1);
        G_n_y(:, 1) = g(x(2:end-1), 0, t_n);
        G_n_y(:, end) = g(x(2:end-1), 1, t_n);
        
        G_n1_y = zeros(m-1, m-1);
        G_n1_y(:, 1) = g(x(2:end-1), 0, t_np1);
        G_n1_y(:, end) = g(x(2:end-1), 1, t_np1);
            
        % solving for intermediate U*
        U_star = M_min_tilde \ ( (M_pl_tilde * U.' + k/(2*h^2) * G_n_y.').' + k/(2*h^2) * G_star);

        % solving for U^{n+1}
        U = (M_min_tilde \ ( (M_pl_tilde * U_star + k/(2*h^2) * G_star).' + k/(2*h^2) * G_n1_y.')).';
    end


    % boundary conditions on final solution 
    U_temp = zeros(m+1);
    U_temp(2:end-1, 2:end-1) = U;
    U = U_temp;
    U(:, 1)   = g(x, 0, t(n));
    U(:, end) = g(x, 1, t(n));
    U(1, :)   = g(0, y, t(n));
    U(end, :) = g(1, y, t(n));
end