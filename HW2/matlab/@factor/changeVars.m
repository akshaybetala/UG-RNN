function f=changeVars(F,fromVar,toVar)
% change the variables over which F is defined
% f=changeVars(F,fromVar,toVar) changes the variables over which F is defined, by
%   swapping each id: fromVar(i) with id: toVar(i)

% (c) Alexander Ihler 2010

  vars=F.v;                % pull out variables for relabeling
  var0=vars;               %   and make a copy for search process
  for i=1:length(fromVar), vars(var0==fromVar(i))=toVar(i); end;
  f=factor(vars,F.t+0);    % may permute indices of the table to preserve ordering property


