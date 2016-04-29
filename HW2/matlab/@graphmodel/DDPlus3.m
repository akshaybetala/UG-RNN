function [logZub,logZlb, tmp] = MAS(fg,order,iBound)
% bound log(Z) via DynaDecompPlus with specified elimination order and i-Bound
% [lnZ+,lnZ-]=MAS(fg,order,iBound) 

errMethod='linf'; decompMethod='L2+hpm'; MAS=false;
randomize = true;
tmp={};

fg=graphmodel(fg.factors); % copy fg to sidestep reference-based updates
Nv=fg.nVar;
epsilon=0; alpha=0;
epsVec = cell(1,length(fg.factors));  for i=1:length(epsVec),epsVec{i}=factor([],0); end;
%fprintf('Current lnZ = %f , err %f\n',log(table(sum(joint(fg),1:16))),epsilon);

step=1; ZERO=factor([],0);
for i=1:length(order)
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
	epsBuc = epsVec(bucketIds); epsVec(bucketIds) = cell(1,length(bucketIds));
  for i=bucketIds, epsVec{i}=factor([],0); end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compute bucket product, splitting too-large factors
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  F=bucket{1}; epsF=epsBuc{1};
  for j=2:Nf, 
    F2=bucket{j}; v2=variables(F2); v=variables(F);
		%v,v2,
    if (length(vunion(v,v2))<=iBound+1)                         % if we're still small enough
      F=F*F2;                                                   %   just include the factor normally
      epsF = epsF+epsBuc{j}; epsBuc{j}=[];			
			%fprintf('Combining in bucket 1...\n'); FP=F*exp(ZERO+epsF),
    else                                                        % otherwise
      Flist = decompProd(F2,{vintersect(v,v2),vdiff(v2,v)},decompMethod);  %   split F2 into smaller components
      F3 = Flist{1}*Flist{2};
      F=F*Flist{1}; [fg,pos]=addFactor(fg,Flist{2});            %   push remainder back
			%fprintf('Splitting; bucket 1:\n'); FP=F*exp(ZERO+epsF),

      tmp{end+1} = log(F2)-log(F3);   % !!! keep track of parameter approximations for test bound
			neweps = abs( log(F2)-log(F3) ) + epsBuc{j}; epsBuc{j}=factor([],0);
			if (~isempty(neweps)) epsVec{pos}=maxmarginal(neweps,diff(v2,v)); end;
			%if (~isempty(neweps)) epsF = epsF + maxmarginal(neweps,v); end;
    end;
  end;
  if (~isempty(epsF))
		FP=sum(F*exp(epsF),X); FM=sum(F*exp(-epsF),X);              % if we have bounds, find eliminated bounds
		F = geomean(FP,FM); epsF = .5*log(FP/FM);                   %   and center our estimate
	else
		F = sum(F,X);                                               % else just eliminate the variable
	end;
%  if (0 && ~isempty(epsF))
%		fprintf('Centering epsF...\n');
%%		epsF,
%		%Fs = sum(F,X); FP1=Fs*exp(max(epsF,X)), FM1=Fs*exp(-max(epsF,X)),
%    q = F/sum(F,X); eta2=sum(epsF*q,X); eta1=max(epsF,X);
%    eta2b=max(epsF,X); eta1b=max(epsF,X);
%		eta2,eta2b,pause;
%		%eta2=maxmarginal(eta2,vars(eta2b));
%    F = F*exp(.5*(eta1-eta2));
%    epsF=.5*(eta1+eta2);
%		FP = sum(F,X) * exp(epsF),
%	%  epsF=max(epsF,X);
%	end;

  %F=sum(F,X);                                                    % eliminate the variable
  [fg,pos]=addFactor(fg,F);                                      % and put the result back into the graph
	if (~isempty(epsF)) epsVec{pos} = max(epsF,X); else epsVec{pos}=[]; end;

end;
      
% Get remaining factors and compute partition function
bucket = getFactor(fg); Nf=length(bucket);          % get all remaining factors
epsF = epsVec{1};
if (Nf>1) F=log(bucket{1}); else F=log(bucket); end;           % get first factor
for j=2:Nf, if (~isempty(bucket{j})),                %  & iterate over the rest
  F=F + log(bucket{j}); epsF=epsF+epsVec{j};
end; end;

%Try to tighten lower bound?
%if (~isempty(epsF)) epsF = marginal( epsF * F / marginal(F,[]) ,[]); end;

logZ=table(logsumexp(F,vars(F)))+alpha;                      % compute normalizer

if (isempty(epsF)) epsF=0; elseif (nvar(epsF)==0) epsF=table(epsF); end;

logZub = logZ+epsF;
logZlb = logZ-epsF;

