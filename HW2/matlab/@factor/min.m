function f = min(F,var)
% min : minimize the factor over a subset of its variables
% f = min(F)      : find the minimum value in F
% f = min(F,vars) : eliminate 'vars' in F by  min_{vars} F

% (c) Alexander Ihler 2010
  if (nargin < 2) f=min(F.t(:)); return; end;

  f=factor; var=uint32(var);
  f.v=vdiff(F.v,var);
  if (isempty(f.v))
    f.t=min(F.t(:));
  else
    mem = vmember(F.v,var); idx=find(mem);
    t=F.t; sz=size(t);
    for i=idx, t=min(t,[],i); end; 
    sz=[sz(~mem) 1]; 
    f.t=reshape(t,sz);
end;

