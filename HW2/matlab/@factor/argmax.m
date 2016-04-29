function [x,val] = argmax(f,which)
% argmax: find a maximizing argument of f(x)
% x = argmax(f,which) finds one or more maximizing arguments of the factor f
%   which = 'first' (fastest), 'all' (all optimizers), or 'random' (a randomly selected one)

% (c) Alexander Ihler 2010

if (nargin<2) which='first'; else which=lower(which); end;
  if (isempty(f.v)), x=[]; val=f.t(1); return; end;   % empty / scalar case

  [val,ind]=max(f.t(:));                  % find the optimal value and first index
  switch(which),                          % depending on passed options,
  case 'first', 
    x=ind2subv(f,ind);                    %   return the first optimum,
  case 'all',
    ind=find(f.t(:)==val);                %   all entries equal to the optimum,
    x=zeros(length(ind),length(f.v));
    for i=1:length(ind), x(i,:)=ind2subv(f,ind(i)); end;
  case 'random',                          
    ind=find(f.t(:)==val);                %   or a randomly chosen entry from that set
    r=fix(length(ind)*rand())+1;
    x=ind2subv(f,ind(r));
  otherwise,                              % anything else is an error
    error(['Unknown option (which=' which ')']);
  end;

