function v = value(gm,tuples)
% return the value of F(x)=\prod fi(x) evaluated at a set of configurations x
% v=value(gm,tuples) : return F(tuple) for each vector of tuples 
%    tuples can be Nt x Nv, where Nv=nvar(gm)

v = ones(size(tuples,1),1);
for f=1:length(gm.factors), if (~isempty(gm.factors{f})),
  v = v .* value(gm.factors{f},tuples(:,vars(gm.factors{f})));
end; end;


