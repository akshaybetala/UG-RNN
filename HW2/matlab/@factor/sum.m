function f = sum(F,var)
% sum out a subset of the variables in F
% f = sum(F)      : compute the sum of all entries in F
% f = sum(F,vars) : eliminate 'vars' in F by  sum_{vars} F(x)

% (c) Alexander Ihler 2010
  if (nargin < 2) f=sum(F.t(:)); return; end;

  f=F; 
  v1=F.v; 
  t=F.t; 
  var=uint32(var);
  v=vdiff(v1,var);
  if (isempty(v))
    t=sum(t(:));
  else
    mem = vmember(v1,var); 
    sz=size(t); 
    for i=find(mem), t=sum(t,i); end;
    sz=[sz(~mem) 1]; 
    t=reshape(t,sz);
  end;
  f.t=t; f.v=v;
