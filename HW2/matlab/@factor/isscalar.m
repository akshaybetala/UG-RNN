function e=isscalar(f)
% isscalar(factor): returns true if the factor is a scalar function (of no variables)

% (c) Alexander Ihler 2010

e = (length(f.t)==1);
