function [logResult,F,origF] = MBE(fg,order,iBound,elimOp,varargin)
% MBE : approximate inference with mini-bucket elimination
% [logResult,belief] = MBE(fg,order,iBound,elimOp [,options]) 
%   order:  a variable elimination order
%   iBound: the clique size not to be exceeded by any mini-bucket; computation is O(d^(iBound+1))
%   elimOp: Determines the elimination operator and bounding direction; one of
%           'sum+' : upper bound the log partition function
%           'sum-' : lower bound the log partition function
%           'max+' : upper bound the maximum probability (MPE or MAP value)
%   logResult:  the log of the approximate result (a bound on the MAP or log partition function)
%   belief:     the un-normalized (sum or max) marginal over any variables not in order
%
%   Options are given as argument pairs: ...name, value...  and include: 
%   'weights', method = 'first', 'distance', 'random'    (for sum+, sum-)
%   'distance', method  -- select heuristic distance type (see factor/distance)
%   'match', true/false -- perform moment matching on buckets before proceeding?
%

switch (elimOp),
  case 'sum+', elim=@sum; elimBound=@max; elimMarg=@marginal;
  case 'sum-', elim=@sum; elimBound=@min; elimMarg=@marginal;
  case 'max+', elim=@max; elimBound=@max; elimMarg=@maxmarginal;
  otherwise, error(['Unknown elimination type ',elimOp]);
end;

%WEIGHTS='equal'; MATCH=1;
WEIGHTS='first';
MATCH=0;
DISTANCE='L1';
for i=1:2:length(varargin),
  switch lower(varargin{i})
  case 'weights', WEIGHTS=varargin{i+1};
  case 'distance', DISTANCE=varargin{i+1};
  case 'match', MATCH=varargin{i+1};
  end;
end;

fg=graphmodel(fg.factors); % copy fg to sidestep reference-based updates
%Nv=max(fg.variables);
Nv=fg.nVar;
lnZ = 0;

original = ones(size(fg.factors));   % which factors are MBE outputs, and which are original?
origF = {};      % save output cliques with original factors

for i=1:length(order)
  %fprintf('===Eliminating %d===\n',order(i));
  X=order(i);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Get bucket for current variable
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  bucketIds = withVariable(fg, X);
  bucket = getFactor(fg, bucketIds );
  Nf = length(bucket); if (Nf==1) bucket={bucket}; end;
  if (Nf==0) continue; end;
  %for j=1:Nf, fg=removeFactor(fg,bucketIds(j)); end;
  fg=removeFactor(fg,bucketIds);

  % keep track of only original buckets also
  origB = bucket; for j=1:Nf, if (original(bucketIds(j))) origB{j}=origB{j}+0; else origB{j}=factor(); end; end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Select allocation into buckets
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  nUsed=Nf; score=zeros(Nf,Nf)-1; tmp=cell(1,Nf);
	if ~isempty(DISTANCE), for j=1:Nf, tmp{j}=elim(bucket{j},X); end; end;
  for pass=1:Nf, 
    if (nUsed < 2) break; end;        % unnecessary for one bucket
    for b1=1:nUsed, for b2=b1+1:nUsed,
      if (length(vunion(variables(bucket{b1}),variables(bucket{b2})))<=iBound)
        if isempty(DISTANCE), score(b1,b2) = 1;  % no distance = scope-based
        else score(b1,b2) = distance( elim( bucket{b1}.*bucket{b2} , X) , tmp{b1}.*tmp{b2} ,DISTANCE);
        end;
      else score(b1,b2)=-1;
      end;
    end; end;
    [mn,ind]=max(score(:)); 
		% !!! find all = mn and choose a random one?
    [b1,b2]=ind2sub(size(score),ind);
    if (mn<0) break; end;             % quit when no more possible
    bucket{b1} = bucket{b1}.*bucket{b2}; tmp{b1}=elim(bucket{b1},X);
    origB{b1} = origB{b1}.*origB{b2}; origB{b2}=origB{nUsed};
    bucket{b2}=bucket{nUsed}; bucket{nUsed}=[]; tmp{b2}=tmp{nUsed}; score(:,nUsed)=-1; nUsed=nUsed-1;
  end;
  for j=1:nUsed, origB{j}=origB{j}.*ones(bucket{j}); end;  % make cliques over full bucket var set
  origF(end+1:end+nUsed) = origB(1:nUsed);

  %if (nUsed>1) fprintf('Approximate; bucket %d not satisfied\n',order(i)); end;
  
  % Do approximate moment matching first?
  if (MATCH)
    var=variables(bucket{1}); for j=2:nUsed, var=vintersect(var,variables(bucket{j})); end;
    for j=1:nUsed, tmp{j}=elimMarg(bucket{j},var); end; F=geomean(tmp{1:nUsed});
    for j=1:nUsed, bucket{j}=bucket{j}.*F./tmp{j}; end;
  end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Weight heuristic:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  switch (WEIGHTS)
  case 'first',   iUseElim = 1;
  case 'distance',                       % use distance heuristic to pick which bucket is sum vs max
    for j=1:nUsed, score(j)=distance( elimBound(bucket{j},X) , elim(bucket{j},X) ,DISTANCE); end;
    [mn,iUseElim]=max(score(1:nUsed));
  case 'random', iUseElim = 1+fix(nUsed*rand);
  end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Eliminate individually within buckets
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for j=1:nUsed, if (~isempty(bucket{j})),
    %fprintf('Mini-factor over: '); fprintf('%d ',uint32(variables(bucket{j}))); fprintf('\n');
    if (j==iUseElim) tmp=elim(bucket{j},X);
    else             tmp=elimBound(bucket{j},X);
    end;
    if (isnan(tmp)||~isfinite(tmp)), X, vars(tmp), pause; end;
	  Ztmp = maxmarginal(tmp,[]); tmp=tmp/Ztmp; lnZ = lnZ+log(Ztmp);
    [fg,pos]=addFactor(fg, tmp);
    original(pos)=0;
    %tmp = sumPower(bucket{j},order(i),invWts(j)); ttmp=table(tmp);
    %if (isnan(tmp)), bucket{j}, tmp, pause; end;
    %if (all(ttmp(:)==0)), 'all zero', bucket{j}, tmp, end;
    %fg=addFactor(fg, sumPower(bucket{j},order(i),invWts(j)) );
  end; end;

end;

% Get remaining factors and compute partition function
bucket = getFactor(fg); Nf=length(bucket);            % get all remaining factors
if (Nf>1) F=bucket{1}; else F=bucket; end;           % get first factor
for j=2:Nf, if (~isempty(bucket{j})),                %  & iterate over the rest
  F=F.*bucket{j};
end; end;
logResult=log(elim(table(F))) + lnZ;                 % compute normalizer

%if (nargout>1) belief=normalize(F); end;                  % compute belief / marginal if desired
