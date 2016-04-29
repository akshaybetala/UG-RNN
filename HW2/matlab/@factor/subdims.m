function d=subdims(f,var)
% subdims : compute the dimensions of a subset of variables
% dim = subdims(f,vars) returns the sizes of variables "vars"
%   only the dimension of those variables that are arguments of f will be returned

  [var,pi]=sort(uint32(var));
  ind=vmember(f.v,var);
  d = size(f.t);
  d = d(ind); d=d(pi);
  if isscalar(d) d=[d,1]; end;
  

