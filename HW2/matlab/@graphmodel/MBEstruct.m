function [logResult,JG] = MBE(fg,order,iBound,elimOp,varargin)
% MBE : approximate inference with mini-bucket elimination
% [logResult,joinGraph] = MBE(fg,order,iBound,elimOp [,options]) 
%   order:  a variable elimination order
%   iBound: the clique size not to be exceeded by any mini-bucket; computation is O(d^(iBound+1))
%   elimOp: Determines the elimination operator and bounding direction; one of
%           'sum+' : upper bound the log partition function
%           'sum-' : lower bound the log partition function
%           'max+' : upper bound the maximum probability (MPE or MAP value)
%   logResult:  the log of the approximate result (a bound on the MAP or log partition function)
%   joinGraph:  a data structure containing the (directed) join graph of the MBE
%
%   Options are given as argument pairs: ...name, value...  and include: 
%   'weights', method = 'first', 'distance', 'random'    (for sum+, sum-) , default first
%   'distance', method  -- select heuristic distance type (see factor/distance), default scope
%   'match', true/false -- perform moment matching on buckets before proceeding? default false
%   'rand', true/false  -- randomize bucket selection process?  default false
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
RANDOM=0;
DISTANCE=[]; % default to scope-based
for i=1:2:length(varargin),
  switch lower(varargin{i})
  case 'weights', WEIGHTS=varargin{i+1};
  case 'distance', DISTANCE=varargin{i+1};
  case 'match', MATCH=varargin{i+1};
  case 'rand', RANDOM=varargin{i+1};
  end;
end;

fg=graphmodel(fg.factors); % copy fg to sidestep reference-based updates
%Nv=max(fg.variables);
Nv=fg.nVar;
lnZ = 0;


JG.cliques={};
JG.cliqElim=[];
JG.eIdx=sparse([]); JG.eSrc=[]; JG.eDst=[];
JG.msg={};
JG.factors={};


Source = -ones(size(fg.factors));   % initially, all factors are "original" => -1

for i=1:length(order)
  %fprintf('===Eliminating %d===\n',order(i));
  X=order(i);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Get bucket for current variable
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  bucketIds = withVariable(fg, X);	% get all factors over X from the temporary model
  bucket = getFactor(fg, bucketIds );
  Nf = length(bucket); if (Nf==1) bucket={bucket}; end;
  if (Nf==0) continue; end; 				% no factors with this variable? should be an error???
  fg=removeFactor(fg,bucketIds);		% remove those factors from the temporary model

  % Keep track of additional information about how the buckets are assembled:
  BucketFactors = cell(1,Nf);  for i=1:Nf, BucketFactors{i} = {bucket{i}}; end;
  BucketSources = cell(1,Nf);  for i=1:Nf, BucketSources{i} = Source(bucketIds(i)); end;
 

  % Do approximate moment matching first?  (Before selecting buckets, eg if distance based)
  if (MATCH && ~isempty(DISTANCE))
    tmp = cell(1,Nf);
    var=variables(bucket{1}); for j=2:Nf, var=vintersect(var,variables(bucket{j})); end;
    for j=1:Nf, tmp{j}=elimMarg(bucket{j},var); end; F=geomean(tmp{1:Nf});  % should use weights?
    for j=1:Nf, bucket{j}=bucket{j} * F / tmp{j}; end;
  end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Select allocation into buckets
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  nUsed=Nf; score=zeros(Nf,Nf)-1; tmp=cell(1,Nf);
	%%%if ~isempty(DISTANCE), for j=1:Nf, tmp{j}=elim(bucket{j},X); end; end;
	if ~isempty(DISTANCE), for j=1:nUsed, tmp{j}=sumPower(bucket{j},X,2); end; end;
  for b1=1:nUsed, for b2=b1+1:nUsed,
    iBoundTmp = max( [iBound+1, nvar(bucket{b1}), nvar(bucket{b2})] );
    %if (iBoundTmp > iBound+1), vars(bucket{b1}),vars(bucket{b2}), pause; end;
    if (length(vunion(variables(bucket{b1}),variables(bucket{b2})))<=iBoundTmp)   
      if isempty(DISTANCE), score(b1,b2) = 1;  % no distance = scope-based
      else score(b1,b2) = distance( elim( bucket{b1}.*bucket{b2} , X) , tmp{b1}.*tmp{b2} ,DISTANCE);
      end;
    else score(b1,b2)=-1;
    end;
  end; end;

	% Join mini-buckets together (by scope or score) until no more joins are possible
  for pass=1:Nf, 
    [mn,ind]=max(score(:));   % find minimum (the first or a random minimum)
    if (RANDOM) ind = find(score(:)==mn); ind=ind(randi(length(ind))); end;
    [b1,b2]=ind2sub(size(score),ind); 				% get which two buckets to merge
    if (b1>b2) btmp=b1; b1=b2; b2=btmp; end;	% (put them in sorted order)
    %fprintf('Got %f : %d %d\n',mn,b1,b2); 
    if (mn<0) break; end;             				% quit when no more merges are possible

    bucket{b1} = bucket{b1}.*bucket{b2}; 			% merge them by taking product
    BucketFactors{b1} = {BucketFactors{b1}{:}, BucketFactors{b2}{:}};  % concatenate list of factors
    BucketSources{b1} = [BucketSources{b1}, BucketSources{b2}];  %  and their "source" info
   

	  %!!! %if ~isempty(DISTANCE), for j=1:nUsed, tmp{j}=sumPower(bucket{j},X,2); end; end;  % prep for re-scoring (???) 
	  if ~isempty(DISTANCE), tmp{b1}=sumPower(bucket{b1},X,2); end;  	% re-elim new bucket for re-scoring
    bucket{b2}=bucket{nUsed}; bucket{nUsed}=[]; tmp{b2}=tmp{nUsed};	% move last bucket to b2's spot and
		score(1:b2,b2)=score(1:b2,nUsed); score(b2,b2+1:nUsed)=score(b2+1:nUsed,nUsed); % copy scores, etc.
    score(:,nUsed)=-1; score(nUsed,:)=-1; score(b2,b2)=-1;
    BucketFactors{b2} = BucketFactors{nUsed}; BucketFactors{nUsed} = {};
    BucketSources{b2} = BucketSources{nUsed}; BucketSources{nUsed} = [];
    nUsed = nUsed - 1;

    for b2=1:nUsed,																% Now, rescore new grouping "b1" with the other groups
      bb1 = min(b1,b2); bb2=max(b1,b2);
      iBoundTmp = max( [iBound+1, nvar(bucket{bb1}), nvar(bucket{bb2})] ); 
      if (length(vunion(variables(bucket{bb1}),variables(bucket{bb2})))<=iBoundTmp)  
        if isempty(DISTANCE), score(bb1,bb2) = 1;  % no distance = scope-based
        else score(bb1,bb2) = distance( elim( bucket{bb1}.*bucket{bb2} , X) , tmp{bb1}.*tmp{bb2} ,DISTANCE);
        end;
      else score(bb1,bb2)=-1;
      end;
    end;
    score(b1,b1)=-1;
  end;

  %Now, add resulting mini-buckets into join-graph structure 
  CliquesAdded = length(JG.cliques) + (1:nUsed);  % adding "nUsed" new mini-bucket cliques (save IDs for later)
  JG.cliqElim(CliquesAdded) = X;  								% these new cliques eliminate "X"
  for j=1:nUsed,
    JG.cliques{end+1} = vars(bucket{j}); 	% clique variables for this mini-bucket
    JG.factors{end+1} = BucketFactors{j}(BucketSources{j}==-1);  % original factors assigned to this clique
    for k=find(BucketSources{j}~=-1), 		% append factors as messages from "src"->"this" 
      thisEdge = length(JG.msg)+1;
      JG.msg{thisEdge} = BucketFactors{j}{k};
      JG.eIdx(BucketSources{j}(k), CliquesAdded(j))=thisEdge;
      JG.eSrc(thisEdge) = BucketSources{j}(k); JG.eDst(thisEdge) = CliquesAdded(j);
      % TODO: repeat for reverse edge with message 1?
    end;
  end;

  %if (nUsed>1) fprintf('Approximate; bucket of var %d not satisfied\n',order(i)); end;
  
  % Do approximate moment matching first? (on buckets, before elimination) (TODO: should use weights!)
  if (MATCH && nUsed > 1)
    var=variables(bucket{1}); for j=2:nUsed, var=vintersect(var,variables(bucket{j})); end;
    for j=1:nUsed, tmp{j}=elimMarg(bucket{j},var); end; F=geomean(tmp{1:nUsed});
    for j=1:nUsed, bucket{j}=bucket{j}.*F / tmp{j}; end;
    % TODO: append messages between mini-buckets (1 to others & others to 1?)
    for j=1:nUsed,  % !!! TODO: HACK just rewrite first factor...
      if (isempty(JG.factors{CliquesAdded(j)})), JG.factors{CliquesAdded(j)} = {factor()}; end;
      JG.factors{CliquesAdded(j)}{1}=JG.factors{CliquesAdded(j)}{1}.*F/tmp{j};
    end;
  end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Weight heuristic:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  switch (WEIGHTS)
  case 'first',   iUseElim = 1;
  case 'distance',                       % use distance heuristic to pick which bucket is sum vs max
    for j=1:nUsed, score(j)=distance( elimBound(bucket{j},X) , elim(bucket{j},X) ,DISTANCE); end;
    [mn,iUseElim]=max(score(1:nUsed));
  case 'random', iUseElim = randi(nUsed); 
  case 'uniform', iUseElim = 1; %randi(nUsed); % neg weights: bucket 1 special?
  end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Eliminate individually within buckets
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %fprintf('Eliminate %d; %d buckets\n',X,nUsed);
  for j=1:nUsed, if (~isempty(bucket{j})),
    %fprintf('Mini-factor over: '); fprintf('%d ',uint32(variables(bucket{j}))); fprintf('\n');
    if (strcmp(WEIGHTS,'uniform'))
      if (strcmp(elimOp,'sum+')) tmp = sumPower(bucket{j},X,nUsed);
      elseif (strcmp(elimOp,'sum-'))
        K=1;
        if (j==iUseElim) tmp = sumPower(bucket{j},X,1/(1+K*(nUsed-1)));%  1*nUsed/(2*nUsed-1));
        else             tmp = sumPower(bucket{j},X,-1/K); 
        end;
      end;
    else
      if (j==iUseElim) tmp=elim(bucket{j},X);
      else             tmp=elimBound(bucket{j},X);
      end;
    end; 
    if (isnan(tmp)||~isfinite(tmp)||max(tmp)==0), X, vars(tmp), bucket{j}, tmp, pause; end;
		Ztmp=max(tmp); tmp=tmp/Ztmp; lnZ = lnZ+log(Ztmp);
    % TODO: doesn't work for heuristic; need "remainder" of logZ 

    [fg,pos]=addFactor(fg, tmp);			% add message into the temporary model
    Source(pos) = CliquesAdded(j);    %   mark them as coming from the newly added clique
  end; end;

end;

% Get remaining factors and compute partition function
bucket = getFactor(fg); Nf=length(bucket);            % get all remaining factors
if (Nf>1) F=bucket{1}; else F=bucket; end;           % get first factor
for j=2:Nf, if (~isempty(bucket{j})),                %  & iterate over the rest
  F=F.*bucket{j};
end; end;
logResult=log(elim(F)) + lnZ;                 % compute normalizer
% "F" returns belief over un-eliminated variables
