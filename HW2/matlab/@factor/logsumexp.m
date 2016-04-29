function f = logsumexp(F,var)
% logsumexp : Log-sum-exponentiate elimination operator
% f = logsumexp(F)      : compute the log sum of table exp(F)
% f = logsumexp(F,vars) : eliminate 'vars' in F by  log[ sum_{vars} exp(F) ]
%

% (c) Alexander Ihler 2010
  if (nargin < 2) m=max(F.t(:)); f=m+log(sum(exp(F.t(:)-m))); return; end;

  f=factor; var=uint32(var);
  f.v=vdiff(F.v,var);
if (1)
  % More numerically stable (normalize by max-marginal)
  Fmax=max(F,var);
  Fmax=substitute(Fmax,-inf,0);
  f = log( sum( exp(F-Fmax) , var) ) + Fmax;
else 
  % Old version (normalize by full table max; slightly faster?)
  Fmax=max(F.t(:)); T=exp(F.t-Fmax);
  if (isempty(f.v))
    f.t=log(sum(T(:)))+Fmax;
  else
    mem = vmember(F.v,var); idx=find(mem); sz=size(T);
    for i=idx, T=sum(T,i); end; 
    %for i=idx(end:-1:1), T=sum(T,i); end; 
    sz=[sz(~mem) 1]; %if (length(sz)<2) sz=[1 sz]; end;
  T=reshape(T,sz);
  f.t=log(T)+Fmax;
  end;
end
