function [alpha_max ] = NS4_compute_max_step(pp, xk, A_ieq, b_ieq )
% Given search direction, compute the maximum step that can be taken in
% that direction
%
% alpha_max satisfies:
% A (x + alpha * pp) <= b
%

%% Parameters

tollerance = 1e-11;
default_max = 10;

%% Computation
% alpha = min (b - A x) ./ A p

% compute terms
b_Ax = b_ieq - A_ieq * xk;
Ap = A_ieq * pp;

assert(all(b_Ax <= 1e-10),'ASSERT FAIL: constraints violated')

% handle finite presision
iToll = abs(Ap) < tollerance;
Ap(iToll) = 0;
iToll = abs(b_Ax) < tollerance;
b_Ax(iToll) = 0;

% discard constraints non-decreasing in p
iNonInc = Ap >= 0;
b_Ax = b_Ax(~iNonInc);
Ap = Ap(~iNonInc);

% compute alpha_max
ratio = b_Ax ./ Ap;
% display(ratio)

ratio(ratio < 0) = [];
if ~isempty(ratio)
    alpha_max = min(ratio);
else
    alpha_max = default_max;
end

% apply upper bound on alpha_max
if alpha_max > default_max
    msg = sprintf('WARNING: alpha_max = %d\nReplaced by %d\n',alpha_max,default_max);
%     disp(msg)
    alpha_max = default_max;
end

end