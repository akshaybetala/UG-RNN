function [gm,pos]=addFactor(gm,f)
% [gm,pos]=addFactor(gm,F) : add factor F to model gm; "pos" is its new factor id in gm
%
  if (gm.nVacant>0)
    pos = gm.vacant(gm.nVacant);
    gm.vacant(gm.nVacant)=0;
    gm.nVacant=gm.nVacant-1;
  else
    %error('No remaining space in factor graph.'); 
    %gm.factors{end+1}=f;   % expand factor list and vacancy queue
    %gm.vacant(end+1)=0;
    pos = length(gm.factors)+1;
  end;

  gm.factors{pos}=f + 0.0;                  % force copy construction
  for j=variables(f),
    gm.vNbrs{j} = vunion(gm.vNbrs{j},pos);
  end;

