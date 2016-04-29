function f = max(F,var)
% max : maximize over a subset of the factor's variablesa
% f = max(F)      : find the maximum value in F
% f = max(F,vars) : eliminate 'vars' in F by  max_{vars} F
%   vars must be in sorted order

% (c) Alexander Ihler 2010
  if (nargin < 2) f=max(F.t(:)); return; end;

  f=F; var=uint32(var);
  f.v=vdiff(F.v,var);
  if (isempty(f.v))
    f.t=max(F.t(:));
  else
    mem = vmember(F.v,var); idx=find(mem);
    sz=size(f.t);
    for i=idx, f.t=max(f.t,[],i); end; 
    sz=[sz(~mem) 1];
    f.t=reshape(f.t,sz);
  end;

