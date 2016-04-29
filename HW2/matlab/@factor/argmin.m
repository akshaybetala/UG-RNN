function [x,val] = argmin(f,which)
% argmin: find a minimizing argument of f(x)
% x = argmin(f,which) finds one or more minimizing arguments of f(x)
%   which = 'first' (fastest), 'all' (all optimizers), or 'random' (a randomly selected one)

% (c) Alexander Ihler 2010

if (nargin<2) which='first'; else which=lower(which); end;

  [val,ind]=min(f.t(:));                  % find optimum value and first index
  switch(which),                          % depending on what the caller requested
  case 'first',                           %   return the first optimum
    x=ind2subv(f,ind);
  case 'all',
    ind=find(f.t(:)==val);                %   a list of all entries equal to the optimum
    x=zeros(length(ind),length(f.v));
    for i=1:length(ind), x(i,:)=ind2subv(f,ind(i)); end;
  case 'random',
    ind=find(f.t(:)==val);                %   or a randomly chosen entry from that set
    r=fix(length(ind)*rand())+1;
    x=ind2subv(f,ind(r));
  otherwise,                              % anything else is an error
    error(['Unknown option (which=' which ')']);
  end;

