function pw = ispairwise(gm)
% ispairwise(gm) : return true if the graph model contains no factors over >2 vars
%

% (c) Alexander Ihler 2010

pw=true;
for i=1:length(gm.factors), if (~isempty(gm.factors{i})),
  if (nvar(gm.factors{i})>2) pw=false; end;
end; end;

