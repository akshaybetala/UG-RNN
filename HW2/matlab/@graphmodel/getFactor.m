function flist=getFactor(gm,i)
% f = getFactor(gm,i) : return a factor or list of factors i
 if (nargin<2), i=1:length(gm.factors); end;  % default: return all factors
 if (isscalar(i)) flist = gm.factors{i};      % scalar index = return a single factor
 else flist = gm.factors(i);                  % vector index = return a cell array of factors
 end;

