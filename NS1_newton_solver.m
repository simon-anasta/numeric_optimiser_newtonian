function [xk ] = NS1_newton_solver(xk, A_ieq, b_ieq, A_eq, Z_new, A_active )
% Run a quasi-newton solver
%

%% Parameters

convergance = 1e-4;

%% Initialize

% clc
format compact

eq_count = 0;

gk = gx(xk);
Fk = Fx(xk);
Bk1 = eye(length(xk) - size(A_active,1) - size(A_eq,1));

%% First Iteration

% initial search direction
[pp, Bk1, Z_new, A_active ] = NS3_compute_search_direction(xk, gk, Bk1, xk, gk, Z_new, A_eq, A_active, -1 );

% compute maximum step length
[alpha_max ] =NS4_compute_max_step(pp, xk, A_ieq, b_ieq );

% conduct line search and locate next iterate
[xk1, Fk1, gk1, alpha ] = NS5_line_search(pp, alpha_max, xk, Fk, gk );

% current slope
slope = sqrt(gk1' * gk1);

% display output
msg = sprintf('      F(x)       |g|       |p|   alpha    a_max  active');
disp(msg)
msg = sprintf(' % 9.1f  % 8.2f  % 8.2f  % 6.2f  % 7.2f  % 3d',Fk1,slope,sqrt(pp'*pp),alpha,alpha_max,size(A_active,1)+size(A_eq,1));
disp(msg)

%% Solver

while slope > convergance
    
    % check for no change in Fk
    if Fk == Fk1
        eq_count = eq_count + 1;
        assert(eq_count < 10,'ASSERT FAIL: to many iterations with no improvement')
    else
        eq_count = 0;
    end
            
    
    % relabel for new iteration
    x_old = xk; xk = xk1;
    g_old = gk; gk = gk1;
    Fk = Fk1;
    Bk = Bk1;
    Z = Z_new;
    
    % update the active set
    [A_active, A_add ] = NS2_update_active_set(xk, A_ieq, b_ieq, alpha, alpha_max, Z, A_eq, A_active );
    
    % compute new search direction
    [pp, Bk1, Z_new, A_active ] = NS3_compute_search_direction(xk, gk, Bk, x_old, g_old, Z, A_eq, A_active, A_add );
    
    % stop if step size is zero or deriviative is too shallow
    if pp' * pp < 1e-12 || norm(Hx(xk,gk)) < 1e-12
        break
    end
    
    % compute maximum step length
    [alpha_max ] = NS4_compute_max_step(pp, xk, A_ieq, b_ieq );
    
    % conduct line search and locate next iterate
    [xk1, Fk1, gk1, alpha ] = NS5_line_search(pp, alpha_max, xk, Fk, gk );
    
    % current slope
    slope = sqrt(gk1' * gk1);
    
    % display output
    msg = sprintf(' % 9.1f  % 8.2f  % 8.2f  % 6.2f  % 7.2f  % 3d',Fk1,slope,sqrt(pp'*pp),alpha,alpha_max,size(A_active,1)+size(A_eq,1));
    disp(msg)
    
end

% display output
msg = sprintf(' % 9.1f  % 8.2f  % 8.2f  % 6.2f  % 7.2f  % 3d',Fk1,slope,sqrt(pp'*pp),alpha,alpha_max,size(A_active,1)+size(A_eq,1));
disp(msg)

end