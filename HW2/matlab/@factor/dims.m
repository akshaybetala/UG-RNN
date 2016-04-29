function dim = dims(f)
% dims: return the variable dimensions (cardiality) for a factor f
% dim = dims(factor) 

% (c) Alexander Ihler 2010

if (isempty(f.v)) dim=[]; return; end;
dim=size(f.t);
if (isscalar(f.v)) dim=dim(1); end;

