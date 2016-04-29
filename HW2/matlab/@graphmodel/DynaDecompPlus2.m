function [logZub,logZlb] = MAS(fg,order,iBound,method)
% bound log(Z) via DynaDecompPlus with specified elimination order and i-Bound
% [lnZ+,lnZ-]=MAS(fg,order,iBound) 

switch(lower(method))
  case 'linf', errMethod='linf'; decompMethod='L2+hpm'; MAS=false;
  %case 'mas',  errMethod = 'mas'; decompMethod='L2+mas'; mas_eps=exp(randn(1)); MAS=true;
  case 'mas',  errMethod = 'mas'; decompMethod='L2+hpm'; mas_eps=exp(randn(1)); MAS=true;
end;
randomize = true;

fg=graphmodel(fg.factors); % copy fg to sidestep reference-based updates
Nv=fg.Nv;
epsilon=0; alpha=0;
epsVec = zeros(1,length(fg.factors));
%fprintf('Current lnZ = %f , err %f\n',log(table(sum(joint(fg),1:16))),epsilon);

for i=1:length(fg.factors), if (~isempty(fg.factors{i})),
  fg.factors{i}=log(fg.factors{i});
  if (MAS), 
    R=table(fg.factors{i}); R=min(R(:));
    if (R<=0) R=R-log(1+mas_eps); fg.factors{i}=fg.factors{i}-R; alpha=alpha+R; end;
  end;
end; end;

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
  F=bucket{1}; epsF=epsVec(bucketIds(1));
  for j=2:Nf, 
    F2=bucket{j}; v2=variables(F2); v=variables(F);
    if (length(vunion(v,v2))<=iBound+1)                         % if we're still small enough
      F=F+F2;                                                   %   just include the factor normally
    else                                                        % otherwise
      if (MAS && table(min(F2,v2))<0)
        %error('Split is less than 1');
        tmp = min(F2,v2) - mas_eps;
        F2=F2-tmp; alpha=alpha+table(tmp);
      end;
      Flist = decompSum(F2,{vintersect(v,v2),vdiff(v2,v)},decompMethod);  %   split F2 into smaller components
      F3 = Flist{1}+Flist{2};
      if (MAS && table(min(F3,v2))<0) 
        %error('Decomp is less than 1'); 
      end;
      switch(errMethod),                                        %   compute error associated with decomposition
        case 'linf', epsilon=epsilon+distance(F2,F3,'Linf');    % calculate L-infinity bound
        case 'mas', epsNew = distance(exp(F3),exp(F2),'mas');              % find the error in this approximation
                    epsNew = (1+epsNew)*(1+epsVec(bucketIds(j)))-1;    % if we already had some error, increase
                    epsilon = max( epsilon, epsNew );           % save overall maximum seen so far
      end;
      F=F+Flist{1}; [fg,pos]=addFactor(fg,Flist{2});            %   push remainder back
      if (MAS) epsF=max(epsF,epsNew); epsVec(pos)=epsNew; end;
    end;
  end;
  F=logsumexp(F,X);                                                    % eliminate the variable
  [fg,pos]=addFactor(fg,F);                                      % and put the result back into the graph
  if (MAS) epsVec(pos)=epsF; end;
  %fprintf('Current lnZ = %f , err %f\n',log(table(sum(joint(fg),1:16)))+alpha,epsilon);
end;
      
% Get remaining factors and compute partition function
bucket = getFactor(fg); Nf=length(bucket);          % get all remaining factors
if (Nf>1) F=bucket{1}; else F=bucket; end;           % get first factor
for j=2:Nf, if (~isempty(bucket{j})),                %  & iterate over the rest
  F=F + bucket{j};
end; end;

logZ=table(logsumexp(F,vars(F)))+alpha;                      % compute normalizer
%logZ=log(sum(table(F)))+alpha;                      % compute normalizer
if (MAS)
  % MAS error measure
  B1 = 1/(1+epsilon) * (logZ+epsilon*alpha);
  B2 = logZ + epsilon*(logZ-alpha);
  logZub = max(B1,B2);
  logZlb = min(B1,B2);
else 
  % Linfinity based method
  logZub = logZ+epsilon; 
  logZlb = logZ-epsilon;
end;

%if (nargout>1) F=normalize(F); end;                  % compute belief / marginal if desired
