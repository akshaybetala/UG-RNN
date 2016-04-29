function f=factor(v,t)
% Factor class constructor function
% f=factor(vars,table)
%   Create a function over discrete variables 'vars', with values specified in 'table'
%   Run 'Contents(factor)' to see a list of accessor and manipulator functions

% (c) Alexander Ihler 2010

if (nargin==0)
  f.v=uint32([]); f.t=1;
elseif (isempty(v))
  f.v=uint32([]); f.t=t+0.0;
else
  [v,ord]=sort(uint32(v));                      % sort copies memory (ensures vars are not shared memory)
  f.v=v;
  if (~isscalar(ord)), f.t=permute(t+0.0,ord);  % ensure table does not share memory also
  else f.t=reshape(t+0.0,numel(t),1);           %  (useful in mex routines; enables writing to reference)
  end;
end;
f=class(f,'factor');

