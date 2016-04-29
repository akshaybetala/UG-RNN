function logZ = MBE(fg,order,iBound)
% logZ=MBE(fg,order,iBound) : upper bound log(Z) via mini-bucket with specified elimination order and i-Bound

fg=graphmodel(fg.factors);  % copy fg to sidestep reference-based updates
Nv=fg.Nv; 
for i=1:Nv-1
  %fprintf('===Eliminating %d===\n',order(i));
  X=order(i);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Get bucket for current variable
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  bucketIds = withVariable(fg, X);
  %fprintf('Getting buckets\n'), bucketIds,
  bucket = getFactor(fg, bucketIds );
  Nf = length(bucket); if (Nf==1) bucket={bucket}; end;
  for j=1:Nf, fg=removeFactor(fg,bucketIds(j)); end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Select allocation into buckets
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for j=2:Nf,
    vars=variables(bucket{j});
    for k=1:j-1
      if (isempty(bucket{k})), bucket{k}=bucket{j}; bucket{j}=[];    % failed, move forward to earliest position
      elseif (length(vunion(variables(bucket{k}),vars))<=iBound+1)      % can add to existing group
        ttmp = table(bucket{k}.*bucket{j});
        %if (all(ttmp==0)) 'All zeros', bucket{k},bucket{j}, pause; end;
        bucket{k}=bucket{k}.*bucket{j}; bucket{j}=[];
      end;
    end;
  end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Set weights to default
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  invWts = zeros(1,Nf)+inf; invWts(1)=1.0;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Eliminate individually within buckets
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for j=1:Nf, if (~isempty(bucket{j})),
  %fprintf('Mini-factor over: '); fprintf('%d ',uint32(variables(bucket{j}))); fprintf('\n');
  [fg,pos]=addFactor(fg, sumPower(bucket{j},order(i),invWts(j)) );
  %fprintf('Putting buckets\n'); pos,
  end; end;

end;

bucket = getFactor(fg);                    % get all remaining factors
Nf = length(bucket); 
if (Nf>1) F=bucket{1}; else F=bucket; end;          % get first factor
for j=2:Nf, if (~isempty(bucket{j})),            % iterate over rest, taking product
  F=F.*bucket{j};
end; end;
logZ=log(sum(table(F)));                  % compute normalizer, log2

