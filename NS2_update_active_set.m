function [A_active, A_add ] = NS2_update_active_set(xk, A_ieq, b_ieq, alpha, alpha_max, ~, A_eq, A_active )
% Check whether or not a new constrain needs to be added to the active set.
% And update the active set if necessary.
%

%% Parameters

tollerance = 1e-11;

%% No constraint needs to be added

if alpha < alpha_max
    % preserve active set
    A_add = -1;
    return
end

%% Constraint needs to be added
% As alpha = alpha_max

A_add = 1;

% identify binding constraints
iBinding = abs(A_ieq * xk - b_ieq) < tollerance;
% at least one more constraint must be binding
assert(sum(iBinding) >= 1 + size(A_active,1),'ASSERT FAIL: Number of binding constraints did not increase')

% list of currently binding constraint
A_ieq = A_ieq(iBinding,:);

% produce list of sorted and new constraints
[A_new, A_active ] = compare_and_sort_constraints(A_active, A_ieq);

while ~isempty(A_new)

    % identify and remove constraints that would not make active set full rank
    indicator = false(size(A_new,1),1);
    for ii = 1:size(A_new,1)
        iRank = rank([A_eq;A_active;A_new(ii,:)]');
        indicator(ii) = iRank == 1 + size(A_active,1) + size(A_eq,1);
    end
    A_new = A_new(indicator,:);    

    % choose new constraint at random (in case there are more than one)
    Ar = rand(size(A_new,1),1);
    iA_new = Ar == max(Ar);
    A_active = [A_active;A_new(iA_new,:)];
    A_new(iA_new,:) = [];
end

end

function [A_new, A_active ] = compare_and_sort_constraints(A_active, A)
% We need to compare every binding constraint with the list of currently
% active constraints, in order to identify the new binding constraint(s).
% To make this easier we keep a sorted order of the active constraints.
%

A_new = zeros(0,size(A,2));

% handle zero active constraints
if isempty(A_active)
    A_new = A;
    return
end

ii = 0;

while ii < size(A,1)
    % increment ii
    ii = ii + 1;
    
    if size(A_active,1) < ii
        % put remaining rows of A aside
        Am = size(A,1) - ii + 1;
        A_new(end+1:end+Am,:) = A(ii:end,:);
        % exit loop
        break
    end
        
    % check for a match
    if all(A(ii,:) == A_active(ii,:))
        % do nothing
    elseif all(A(ii,:) == A_active(end,:))
        % reorder A_active
        row = A_active(end,:);
        A_active(end,:) = [];
        A_active = [A_active(1:(ii-1),:);row;A_active(ii:end,:)];
    else
        % put row ii of A aside
        A_new(end+1,:) = A(ii,:);
        % remove row ii of A
        A(ii,:) = [];
        % step ii back one step to check this row again
        ii = ii - 1;
    end
end

end