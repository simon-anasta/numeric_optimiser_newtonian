function [xk1, Fk1, gk1, alpha ] = NS5_line_search(pp, alpha_max, xk, Fk, gk )
% Does line search in direction pp for new iteation:
% xk1 = xk + alpha * pp
%

%% Parameters

alpha_min1 = 0.1 / sqrt(pp' * pp);
alpha_min2 = min(0.001, 0.001 / sqrt(pp' * pp));
rhi_max = 0.2;
max_iterations = 20;
tollerance = 1e-11;

% assert descent direction
assert(gk' * pp <= tollerance,'ASSERT FAIL: p is not a descent direction')

%% Check alpha_min

% specify intended alpha_min
alpha_min = min(alpha_min1, alpha_min1 * alpha_max);
% compute function at step size
xk_test = xk + alpha_min * pp;
Fk_test = Fx(xk_test);
gk_test = gx(xk_test);
% check for acceptability
% if gk_test' * pp <= 0 && Fk_test < Fk then everything is fine
if gk_test' * pp > 0 || Fk_test > Fk
    % reduce alpha_min
    alpha_min = min(alpha_min2, alpha_min2 * alpha_max);
end

% display(alpha_min)

%% Aim for step size 1

if alpha_max >= 1
    % compute function at step size
    xk_test = xk + pp;
    Fk_test = Fx(xk_test);
    gk_test = gx(xk_test);
    % check gradient is small enough
    if abs(gk_test' * pp) <= - rhi_max * gk' * pp && Fk_test < Fk
        % accept step size and return solution
        alpha = 1;
        xk1 = xk_test;
        gk1 = gk_test;
        Fk1 = Fk_test;
        return
    end
end

%% Check step size alpha_max

% computer function at step size
xk_test = xk + alpha_max * pp;
Fk_test = Fx(xk_test);
gk_test = gx(xk_test);
% check gradient is negative
if gk_test' * pp < 0 && Fk_test <= Fk
    % accept step size and return solution
    alpha = alpha_max;
    xk1 = xk_test;
    gk1 = gk_test;
    Fk1 = Fk_test;
    return
end

%% Use interval reduction to search for appropriate alpha

interval = [alpha_min,alpha_max];
shallow = false;
alpha_test = mean(interval);
counter = 0;

while ~shallow
    
%     display(interval)
    counter = counter + 1;
    
    % computer function at step size
    xk_test = xk + alpha_test * pp;
    Fk_test = Fx(xk_test);
    gk_test = gx(xk_test);
    if abs(gk_test' * pp) <= - rhi_max * gk' * pp && Fk_test < Fk
        % accept step size and return solution
        shallow = true;
    elseif Fk_test > Fk
        % reduce alpha
        interval(2) = alpha_test;
        if interval(2) > 1
            alpha_test = sqrt(prod(interval));
        else
            alpha_test = mean(interval);
        end
        
    elseif Fk_test <= Fk && gk_test' * pp < 0
        % increase alpha
        interval(1) = alpha_test;
        alpha_test = mean(interval);
        
    elseif Fk_test <= Fk && gk_test' * pp > 0
        % reduce alpha
        interval(2) = alpha_test;
        if interval(2) > 1
            alpha_test = sqrt(prod(interval));
        else
            alpha_test = mean(interval);
        end
        
    end
    
    % check for no solution
    if range(interval) < 1e-8
        assert(false,'ASSERT FAIL: interval is narrow without convergence')
    end
    
    % check for too many iterations
    if counter > max_iterations
        assert(false,'ASSERT FAIL: step size not found after %d iterations',max_iterations)
    end
end

alpha = alpha_test;
xk1 = xk_test;
gk1 = gk_test;
Fk1 = Fk_test;

end