function M = matrix(F,var1,var2)
% matrix : convert factor to a two-dimensional "pairwise" representation for two variable sets
% M = matrix(F,var1 [,var2]) converts a multivariate factor into a two-dimensional matrix
%   M is a |v1|x|v2| matrix of factor values
%   vunion(var1,var2) must equal vars(F), but can be overlapping sets
%   if var2 is not passed, it is assumed to be vars(F)\var1   (and is slightly faster)

% (c) Alexander Ihler 2010

if (nargin < 3)
  % if we have only two arguments, var1 & var2 are disjoint (var2=F.v - var1) and things are easier:
  var1 = sort(uint32(var1));
  ind = vmember(F.v,var1);
  M = F.t; sz=size(M); D1=prod(sz(ind)); D2=prod(sz(~ind));
  if (length(ind)>1) M=permute(M, [find(ind) find(~ind)]); end;
  M = reshape(M , D1,D2);
  %M = reshape(permute(M, [find(ind) find(~ind)]) , D1,D2);
else
  % if we might be requesting overlapping variable sets
  %  error('Overlapping version not yet implemented');
  vars = F.v; mx=max(vars); dim = dims(F);
  var1 = sort(uint32(var1)); var2 = sort(uint32(var2)); % var12=vintersect(var1,var2);
  mem1 = vmember(vars,var1); mem2 = vmember(vars,var2); mem12 = mem1 & mem2;
  Fid=factor();
  for i=find(mem12), Fid=Fid * factor([vars(i),mx+i],eye(dim(i))); end;
  Fid = F * Fid;
  M = Fid.t; D1=prod(dim(mem1)); D2=prod(dim(mem2));
  Nd=length(dim); idx=1:Nd; idx(mem12)=Nd+1:Nd+sum(mem12);
  M = reshape( permute(M, [find(mem1) idx(mem2)]) ,D1,D2);
end;

