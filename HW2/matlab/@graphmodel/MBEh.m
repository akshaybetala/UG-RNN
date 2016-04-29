function [logZ,F] = MBE(fg,order,iBound)
% logZ=MBE(fg,order,iBound) : upper bound log(Z) via mini-bucket with specified elimination order and i-Bound

%WEIGHTS='equal'; MATCH=1;
WEIGHTS='first';
MATCH=0;

fg=graphmodel(fg.factors); % copy fg to sidestep reference-based updates
Nv=fg.nVar;
for i=1:Nv-1
  %fprintf('===Eliminating %d===\n',order(i));
  X=order(i);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Get bucket for current variable
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  bucketIds = withVariable(fg, X);
  bucket = getFactor(fg, bucketIds );
  Nf = length(bucket); if (Nf==1) bucket={bucket}; end;
  for j=1:Nf, fg=removeFactor(fg,bucketIds(j)); end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Select allocation into buckets
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  nUsed=Nf; score=zeros(Nf,Nf)-1; tmp=cell(1,Nf);
  for j=1:Nf, tmp{j}=sum(bucket{j},X); end;
  for pass=1:Nf, 
    for b1=1:nUsed, for b2=b1+1:nUsed,
      if (length(vunion(variables(bucket{b1}),variables(bucket{b2})))<=iBound)
        score(b1,b2) = distance( sum( bucket{b1}.*bucket{b2} , X) , tmp{b1}.*tmp{b2} ,'L1');
      else score(b1,b2)=-1;
      end;
    end; end;
    [mn,ind]=max(score(:)); [b1,b2]=ind2sub(size(score),ind); %score, [b1,b2],
    if (mn<0) break; end;             % quit when no more possible
    bucket{b1} = bucket{b1}.*bucket{b2}; tmp{b1}=sum(bucket{b1},X);
    bucket{b2}=bucket{nUsed}; bucket{nUsed}=[]; tmp{b2}=tmp{nUsed}; score(:,nUsed)=-1; nUsed=nUsed-1;
  end;

  %if (nUsed>1) fprintf('Approximate; bucket %d not satisfied\n',order(i)); end;
  
  % Do approximate moment matching first?
  if (MATCH)
    var=variables(bucket{1}); for j=2:nUsed, var=vintersect(var,variables(bucket{j})); end;
    for j=1:nUsed, tmp{j}=marginal(bucket{j},var); end; F=geomean(tmp{1:nUsed});
    for j=1:nUsed, bucket{j}=bucket{j}.*F./tmp{j}; end;
  end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Weight heuristic:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  switch (WEIGHTS)
  case 'first', invWts = zeros(1,nUsed)+inf; invWts(1)=1.0;  % first bucket gets 1, others get max
  case 'equal', invWts = ones(1,nUsed).*nUsed;        % all buckets get equal weight
  case 'bestmax',                       % find the greedily "best" bucket to get weight 1
    for j=1:nUsed, score(j)=distance( max(bucket{j},X) , sum(bucket{j},X) ,'L1'); end;
    [mn,ind]=max(score(1:nUsed));
    invWts = zeros(1,Nf)+inf; invWts(ind)=1.0;
  case 'random', invWts=zeros(1,nUsed)+inf; invWts(1+fix(nUsed*rand))=1;
  end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Eliminate individually within buckets
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for j=1:Nf, if (~isempty(bucket{j})),
    %fprintf('Mini-factor over: '); fprintf('%d ',uint32(variables(bucket{j}))); fprintf('\n');
    tmp = sumPower(bucket{j},order(i),invWts(j)); ttmp=table(tmp);
    %if (isnan(tmp)), bucket{j}, tmp, pause; end;
    %if (all(ttmp(:)==0)), 'all zero', bucket{j}, tmp, end;
    fg=addFactor(fg, sumPower(bucket{j},order(i),invWts(j)) );
  end; end;

end;

% Get remaining factors and compute partition function
bucket = getFactor(fg); Nf=length(bucket);            % get all remaining factors

%for j=1:Nf, vars(bucket{j}), end; pause;

if (Nf>1) F=bucket{1}; else F=bucket; end;           % get first factor
for j=2:Nf, if (~isempty(bucket{j})),                %  & iterate over the rest
  F=F.*bucket{j};
end; end;
logZ=log(sum(table(F)));                      % compute normalizer

if (nargout>1) F=normalize(F); end;                  % compute belief / marginal if desired
