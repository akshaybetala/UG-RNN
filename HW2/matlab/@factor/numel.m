function n = numel(F)
% return the number of elements in the specification of the factor's values
% n = numel(f) : returns the number of elements in the table for f (the number of configurations of variables)

% (c) Alexander Ihler 2010

n = prod(size(F.t));
