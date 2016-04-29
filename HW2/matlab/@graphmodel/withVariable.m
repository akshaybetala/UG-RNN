function nbrs = withVariable(gm,list)
% nbrs = withVariable(fg,list) : return the set of factors that include any variables in 'list'
  nbrs = uint32([]);
  for j=1:length(list),
    nbrs=vunion(nbrs,gm.vNbrs{list(j)});
  end;
