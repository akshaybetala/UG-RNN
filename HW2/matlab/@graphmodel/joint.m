function [F,logZ]=joint(gm)
% [F,logZ] = joint(graphmodel) : compute the full joint distribution F of the graphical model
%  optionally return the log partition function as well

F=gm.factors{1};
for i=2:length(gm.factors), if (~isempty(gm.factors{i})), F=F.*gm.factors{i}; end; end;
if (nargout>1) logZ=log(table(sum(F,variables(F)))); end;

