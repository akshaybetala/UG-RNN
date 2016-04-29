function f=condition(F,vars,vals)
% compute the conditional function F(X|Xv=V)
% f=condition(F,vars,vals) : compute the conditional function f=F(x|x_vars=vals)
%  vars must be in sorted order

% (c) Alexander Ihler 2010

  idx=cell(1,length(F.v)); vars=uint32(vars);      % set up input and cell array vars
  cond=vmember(F.v,vars);                          % get the conditioned indices
  vals=vals(vmember(vars,F.v));                    % make sure there aren't any extra variables passed
  idx(cond)=num2cell(vals); idx(~cond)={':'};      % make a cell array of indices, e.g. (v1,:,:,v2,...)
  f=factor(vdiff(F.v,vars),squeeze(F.t(idx{:})));  % make new, conditional factor
