function [logZub,logZlb] = MAS(fg,order,iBound,method)
% bound log(Z) via DynaDecompPlus with specified elimination order and i-Bound
% [lnZ+,lnZ-]=MAS(fg,order,iBound) 

switch(lower(method))
  case 'linf', errMethod='linf'; decompMethod='L2+hpm'; MAS=false;
  case 'mas',  errMethod = 'mas'; decompMethod='L2+hpm'; mas_eps=exp(randn(1)); MAS=true;
end;
randomize = true;

fg=graphmodel(fg.factors); % copy fg to sidestep reference-based updates
Nv=fg.nVar;
nF=length(fg.factors);
epsilon=0; alpha=0;
epsVec = cell(1,nF); for i=1:nF, epsVec{i}=factor([],0); end;
%fprintf('Current lnZ = %f , err %f\n',log(table(sum(joint(fg),1:16))),epsilon);

for i=1:Nv-1
  %fprintf('===Eliminating %d===\n',order(i));
  X=order(i);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Get bucket for current variable
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  bucketIds = withVariable(fg, X);
  if (length(bucketIds)==0) continue; end;
  bucket = getFactor(fg, bucketIds );
  Nf = length(bucket); if (Nf==1) bucket={bucket}; end;
  for j=1:Nf, fg=removeFactor(fg,bucketIds(j)); end;
  if (randomize) pi=randperm(Nf); bucket=bucket(pi); bucketIds=bucketIds(pi); end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compute bucket product, splitting too-large factors
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  F=bucket{1}; epsF=epsVec{bucketIds(1)}; epsVec{bucketIds(1)}=factor([],0);
  for j=2:Nf, 
    F2=bucket{j}; v2=variables(F2); v=variables(F);
    if (length(vunion(v,v2))<=iBound)                           % if we're still small enough
      F=F*F2;                                                   %   just include the factor normally
			epsF=epsF+epsVec{bucketIds(j)}; epsVec{bucketIds(j)}=factor([],0);
    else                                                        % otherwise
      Flist = decompProd(F2,{vintersect(v,v2),vdiff(v2,v)},decompMethod);  %   split F2 into smaller components
      F3 = Flist{1}*Flist{2};
      err = abs(log(F2)-log(F3)) + epsVec{bucketIds(j)}; epsVec{bucketIds(j)}=factor([],0);
      F=F*Flist{1}; [fg,pos]=addFactor(fg,Flist{2});            %   push remainder back
			epsVec{pos} = maxmarginal(err,vars(Flist{2}));
    end;
  end;
  F=sum(F,X);                                                    % eliminate the variable
  [fg,pos]=addFactor(fg,F);                                      % and put the result back into the graph
	epsVec{pos}=max(epsF,X);
  %fprintf('Current lnZ = %f , err %f\n',log(table(sum(joint(fg),1:16)))+alpha,epsilon);
end;
      
% Get remaining factors and compute partition function
bucket = getFactor(fg); Nf=length(bucket);          % get all remaining factors
if (Nf>1) F=log(bucket{1}); else F=log(bucket); end;           % get first factor
eps = epsVec{1};
for j=2:Nf, if (~isempty(bucket{j})),                %  & iterate over the rest
  F=F + log(bucket{j});
	eps = eps + epsVec{j};
end; end;

logZ=table(logsumexp(F,vars(F)))+alpha;                      % compute normalizer
%logZ=log(sum(table(F)))+alpha;                      % compute normalizer
epsilon = table(max(eps,vars(eps)));

logZub = logZ+epsilon; 
logZlb = logZ-epsilon;

%if (nargout>1) F=normalize(F); end;                  % compute belief / marginal if desired
