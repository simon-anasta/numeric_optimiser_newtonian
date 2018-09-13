function [pp, Bk1, Z, A_active_new ] = NS3_compute_search_direction(xk, gk, ~, ~, ~, Z, A_eq, A_active, A_add )
% Given active set checks whether or not to move off constraints,
% identifies constraints to move off and computes search direction.
%

%% Parameters

lambda_gamma = 0.1;
tollerance = 1e-12;

%% Initial computing

% If moved onto a new constraint
if A_add == 1
    % Recompute Z
    [q,~] = qr([A_eq;A_active]');
    deg_freedom = length(xk) - size(A_active,1) - size(A_eq,1);
    Z = q(:,end-deg_freedom+1:end);
    
    % Check reduced basis Z is orthogonal to active constraints
    assert(all(all(abs([A_eq;A_active] * Z) <= tollerance)),'ASSERT FAIL: reduced basis not orthogonal to constraints (1)')
end

% Compute Hessian matrix
Bk1 = Hx(xk,gk);

%% Compute lambda_max

A_active_new = A_active;

if size(A_eq,1) + size(A_active_new,1) == size(xk,1)
    lambda_max = 0;
else
    lambda_max = -lambda_gamma * norm(Z' * gk);
    if abs(lambda_max) < tollerance
        lambda_max = 0;
    end
end

%% Identify whether or not to move off constraints

% compute search direction
[qq, dd ] = search_direction(gk, Z, Bk1 );
pp = qq + dd;

% If there are inequality constraints in the active set
% And we did not move onto a constraint last iteration
if ~isempty(A_active_new) && A_add ~= 1
    % find lambda that minimizes |A' lambda - (g + H p)|^2
    lambda = [A_eq;A_active_new]' \ (gk + Bk1 * qq);
    lambda(abs(lambda) < tollerance) = 0;
    % keep only lambdas corresponding to active set
    lambda(1:size(A_eq,1)) = [];
    
    % select constraint to move off
    iLambda = lambda <= lambda_max & lambda == min(lambda) & min(lambda) < 0;
    
    % prevent moving off more than one constraint
    if sum(iLambda) > 1
        lambda(iLambda) = lambda(iLambda) - rand(sum(iLambda),1);
        iLambda = lambda <= lambda_max & lambda == min(lambda) & min(lambda) < 0;
    end
    
    % compute new A and Z
    A_active_new(iLambda,:) = [];
    z = null([A_eq' A_active_new' Z]');
    Z = [Z z];
    assert(size(z,2) <= 1,'ASSERT FAIL: null space increases by extra dimensions')
    
    % reduced basis Z must be orthogonal to constraints
    assert(all(all(abs([A_eq;A_active_new] * Z) <= tollerance)),'ASSERT FAIL: reduced basis not orthogonal to constraints (2)')

    % return early if moving off no constraints
    if sum(iLambda) == 0
        return
    end
    
    % recompute search direction
    [qq, dd ] = search_direction(gk, Z, Bk1 );
    pp = qq + dd;
end

end

%% Sub function : Modified Cholesky Factorization

function [L, E ] = mod_cholesky(G)
% Given a symmetric matrix, carry out modified Cholesky factorization:
% - find E such that |E| is small and G + E is positive definite
% - and find lower tirangular L such that L*L' = G+E
%
% Inspired by mchol1 from Matlab File Exchange
%  Author: Brian Borchers (borchers@nmt.edu)
%  Reference: Gill, Murray, and Wright, "Practical Optimization", p111.
%

%% Parameters

n = size(G,1);

%% Optimal Beta

% off diagonals
zi = max(max(abs(G-diag(diag(G)))));
% diagonals
gamma = max(diag(G));
% tollerance
tollerance = 1e-12;
% Beta squared
nu = max([1,sqrt(n^2-1)]);
beta = sqrt(max([gamma * 1.01, zi/nu, tollerance]));

%% Modified Cholesky

A = G;
R = zeros(size(A));
E = zeros(n,1);

for kk = 1:n
    % ensure diagonal element is large enough
    mu_k = max(abs(A(kk,kk+1:n)));
    R(kk,kk) = max([tollerance, sqrt(abs(A(kk,kk))), mu_k / beta]);
    E(kk) = R(kk,kk)^2 - A(kk,kk);
    % compute off-diagonal entries
    R(kk,kk+1:n) = A(kk,kk+1:n) ./ R(kk,kk);
    for jj = kk+1:n
        A(kk+1:jj,jj) = A(kk+1:jj,jj) - R(kk,jj) .* R(kk,kk+1:jj)';
    end
end

%% Prepare output

E = diag(E);
L = R';

end

%% Sub function : Compute Step Direction

function [qq, dd ] = search_direction(gk, Z_new, Bk1 )
% Given the active set and current point compute the descent direction and
% any direction of negative curvature.
%

%% Newton search direction

% Reduced Hessian
Bk1 = Z_new' * Bk1 * Z_new;

% pp solves Z'BZ pp = -Z'g
gz = -Z_new' * gk;
if isempty(Bk1) || norm(Bk1) < 1e-8
    q_z = gz;
else
    L = mod_cholesky(Bk1);
    tmp = L \ gz;
    q_z = L' \ tmp;
end

% pp = Z p_z
qq = Z_new * q_z;

assert(gk' * qq <= 1e-10,'ASSERT FAIL: Newton direction is not a descent direction')

%% Direction of negative curvature

if isempty(Bk1) || norm(Bk1) < 1e-10
    % no direction of negative curvature
    dd = zeros(size(qq));
else
    % compute eigen values
    [evectors,evalues] = eig(Bk1);
    evalues = diag(evalues);
    
    if min(evalues) >= 0
        % no direction of negative curvature
        dd = zeros(size(qq));
    else
        % select eigenvector corresponding to smallest eigenvalue
        iMin = evalues == min(evalues);
        d_z = evectors(:,iMin);
        % in case multiple minimum eigenvalues
        d_z = d_z(:,1);
        % nomalize
        d_z = d_z * abs(min(evalues));
        % convert to full space
        dd = Z_new * d_z;
        
        % check whether we want dd or -dd
        if gk' * dd > 1e-10
            dd = -dd;
        end
        assert(gk' * dd <= 1e-10,'ASSERT FAIL: direction of negative curvature is not a descent direction')
    end
end


end

