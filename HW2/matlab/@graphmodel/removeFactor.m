function gm=removeFactor(gm,ilist)

for i=uint32(ilist),
  %if (isempty(fg.factors{i})) return; end;    % if already empty, do nothing
  %if (gm.factors{i}==1) return; end;

  for j=variables(gm.factors{i}),              % otherwise, remove from adjacency structure
    gm.vNbrs{j} = vdiff(gm.vNbrs{j},i);
  end;
	%modCell(gm.factors,i,factor());
  gm.factors{i}=factor();                      % clear entry
  
  gm.nVacant = gm.nVacant+1;                  % add to list of available spots
  gm.vacant(gm.nVacant)=i;
  %fg.fEmpty=vunion(fg.fEmpty,i);          % and add to list of available spots
end;
  
