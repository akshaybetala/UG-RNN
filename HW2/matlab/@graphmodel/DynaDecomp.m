function [logZub,logZlb] = MAS(fg,order,iBound,method)
% upper bound log(Z) via DynaDecomp with specified elimination order and i-Bound
% [lnZ+,lnZ-]=MAS(fg,order,iBound) 

switch(lower(method))
  case 'linf', errMethod='linf'; decompMethod='L2+hpm'; MAS=false;
  %case 'mas',  errMethod = 'mas'; decompMethod='L2+mas'; mas_eps=exp(randn(1)); MAS=true;
  case 'mas',  errMethod = 'mas'; decompMethod='L2+hpm'; mas_eps=exp(randn(1)); MAS=true;
end;
randomize = true;

fg=graphmodel(fg.factors); % copy fg to sidestep reference-based updates
nVar=fg.nVar;
epsilon=0; alpha=0;
epsVec = zeros(1,length(fg.factors));
%fprintf('Current lnZ = %f , err %f\n',log(table(sum(joint(fg),1:16))),epsilon);

if (MAS), for i=1:length(fg.factors), if (~isempty(fg.factors{i})),    % if MAS error methods,
  R = table(min(fg.factors{i},vars(fg.factors{i})));              % make sure all factors are > 1
  if (R<=1) R=R/exp(1+mas_eps); fg.factors{i}=fg.factors{i}/R; alpha=alpha+log(R); end;
end; end; end;

for i=1:nVar-1
  %fprintf('===Eliminating %d===\n',order(i));
  X=order(i);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Get bucket for current variable
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  bucketIds = withVariable(fg, X);
  bucket = getFactor(fg, bucketIds );
  Nf = length(bucket); if (Nf==1) bucket={bucket}; end;
  for j=1:Nf, fg=removeFactor(fg,bucketIds(j)); end;
  if (randomize) pi=randperm(Nf); bucket=bucket(pi); bucketIds=bucketIds(pi); end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Eliminate X and return subject to size constraints 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  F=bucket{1};
  for j=2:Nf, F=F*bucket{j}; end;
  F=sum(F,X);                % eliminate X and check the size
  nvF = nvar(F);
  if (nvF <= iBound+1)
    fg = addFactor(fg,F);    % if it's not too big just add it back
  else                      % otherwise:
    v=variables(F); pi = v(randperm(nvF));                             % choose random subsets of variables
    nNew = ceil(nvF/(iBound+1)); nPer = (nvF/nNew)-1; sets=cell(1,nNew);  % identify smaller sets
    jj=1; for j=1:nNew, sets{j}=sort(pi( round(jj):min(nvF,round(jj+nPer)) )); jj=round(jj+nPer+1); end;
    decomp = decompProd(F, sets, decompMethod);                         % split F into these subsets
    F2=decomp{1}; for j=2:length(decomp), F2=F2*decomp{j}; end;        % measure error from split
    if (MAS && table(min(F2,v))<1) error('Decomp is less than 1'); end;
    switch(errMethod),
      case 'linf', epsilon = epsilon + distance(log(F),log(F2),'Linf');% "safe" calculation in case not centered
      case 'mas', epsNew = distance(F2,F,'mas');                       % find the error in this approximation
                  epsNew = (1+epsNew)*prod(1+epsVec(bucketIds))-1;     %   if we already had some error, increase
                  epsilon = max( epsilon, epsNew );                   % save overall maximum seen so far
    end;
    for j=1:length(decomp), 
      [fg,pos]=addFactor(fg, decomp{j});          % push back the remaining factors
      if (MAS) epsVec(pos) = epsNew; end;
    end;
  end;
  %fprintf('Current lnZ = %f , err %f\n',log(table(sum(joint(fg),1:16)))+alpha,epsilon);
end;

% Get remaining factors and compute partition function
bucket = getFactor(fg); Nf=length(bucket);          % get all remaining factors
if (Nf>1) F=log(bucket{1}); else F=log(bucket); end;           % get first factor
for j=2:Nf, if (~isempty(bucket{j})),                %  & iterate over the rest
  F=F + log(bucket{j});
end; end;

logZ=table(logsumexp(F,vars(F)))+alpha;                      % compute normalizer
%logZ=log(sum(table(F)))+alpha;                      % compute normalizer
if (MAS)
  % MAS error measure
  B1 = 1/(1+epsilon) * (logZ+epsilon*alpha);
  B2 = logZ+epsilon*(logZ-alpha);
  logZub = max(B1,B2);
  logZlb = min(B1,B2);
else
  % Linfinity based method
  logZub = logZ+epsilon;
  logZlb = logZ-epsilon;
end;

%if (nargout>1) F=normalize(F); end;                  % compute belief / marginal if desired
