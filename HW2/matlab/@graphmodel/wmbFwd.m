function [gm, lnZ] = wmbFwd(gm, updateThetas, updateWeights)
% [gm, lnBound]=wmbFwd(gm [,options]) : compute forward messages for weighted mini-bucket elimination

if (~strcmp(gm.Alg.name,'WMB'))
  fprintf('Graphical model has not been specialized to the WMB algorithm\n');
  return;
end;

lnZ = 0.0;
if (nargin < 3) updateWeights=true; end;  % default settings
if (nargin < 2) updateThetas=true;  end;

if (islogical(updateWeights)) updateWeights = 0.5; end;
if (islogical(updateThetas)) updateWeights = 0.1; end;

for i=1:length(gm.Alg.order),  % for each bucket (= collection of minibuckets)
  X = gm.Alg.order(i); %minibucket(i).var;     % variable to be eliminated
  nodes = gm.Alg.nodes{i}; 
  nNodes = length(gm.Alg.nodes{i});  % # of minibucket partitions in bucket
  bel = cell(1,nNodes);             % allocate storage for beliefs

  if (length(gm.Alg.nodes{i})>1) % if more than one mini-bucket partition:

    if (updateThetas || updateWeights)                       % if match, compute beliefs
      for j=1:nNodes,
        n = gm.Alg.nodes{i}(j);                   % compute belief at node n:
        bel{j}=gm.Alg.theta{n}+gm.Alg.msgBwd{n}; % theta+bwd+\sum fwd
        for c=gm.Alg.children{n}, bel{j}=bel{j}+gm.Alg.msgFwd{c}; end;
        bel{j}=bel{j}*(1.0/gm.Alg.wt(n));              % and power (TODO?)
        bel{j}=bel{j} - logsumexp(bel{j});                   % and normalize
      end;
    end;

    if(updateThetas)
      % Update theta (parameter allocation)
      % (TODO): save factor marginals somewhere here...
      for m=1:length(gm.Alg.match{i}), 
        match = gm.Alg.match{i}{m};
        wTot=0; 
        for c=match, wTot=wTot+gm.Alg.wt(c); end;
        vAll=gm.Alg.clique{match(1)};
        for c=match, vAll=vintersect(vAll,gm.Alg.clique{c}); end;
        delta=cell(1,nNodes); bavg=log(factor());
        for c=match,
          j = find( gm.Alg.nodes{i} == c );
          delta{j} = logsumexp(bel{j}, vdiff(vars(bel{j}),vAll));  % get marginal over match vars
          bavg = bavg + delta{j}*(gm.Alg.wt(c) / wTot); % & weighted avg
          %OR: compute actual average (for proj gradient)  (TODO)
        end;
        damp = updateThetas; %.5; % TODO
        for c=match,
          j = find( gm.Alg.nodes{i} == c );
          delta{j} = (bavg - delta{j});                        % compute update & apply to
          bel{j} = bel{j} + delta{j}*damp;                     %   belief (TODO?)
          bel{j}=bel{j} - logsumexp(bel{j});                   %   and normalize (TODO?)
          gm.Alg.theta{c} = gm.Alg.theta{c} + delta{j}*gm.Alg.wt(c)*damp; % and theta
        end;
      end; % (list of matches)
    end;   % (update thetas)

    allWeights = [gm.Alg.wt(nodes)];
    isMax = sum(allWeights)<.9;        % doesn't sum to 1 => max+
    isLower = any(allWeights<0);       % negative weights => sum-
    if(updateWeights && ~isMax)
      wwstep = updateWeights; %.1;  % (TODO) step size for weights
      Havg = 0; nPos=0; H=cell(1,nNodes); wTot=0;
      if (any(allWeights<0) && sum(allWeights>0)~=1) error('Weights set incorrectly?'); end;
      for j=1:nNodes,
        n = gm.Alg.nodes{i}(j);       % compute conditional entropy at node n:
        H{j} = -sum( exp(bel{j}) * (bel{j}-logsumexp(bel{j},X)) );
        if (~isLower)
          Havg = Havg + gm.Alg.wt(n)*H{j};   % upper bound: track average entropy
        elseif (gm.Alg.wt(n) > 0)            % lower: track only positive weight entry
          Havg = H{j}; nPos=n; 
        end;  
      end;

      for j=1:nNodes,
        n = gm.Alg.nodes{i}(j);       % take a step in the gradient direction 
        if (~isLower)
          gm.Alg.wt(n) = gm.Alg.wt(n) * exp(-wwstep*gm.Alg.wt(n)*(H{j}-Havg));
          wTot = wTot + gm.Alg.wt(n);      % and compute the weights of the new point
        elseif (gm.Alg.wt(n) < 0), 
          gm.Alg.wt(n) = gm.Alg.wt(n) * exp(wwstep*gm.Alg.wt(n)*(H{j}-Havg));
          wTot = wTot + gm.Alg.wt(n);
        end;
      end;                                       
      if (~isLower),                             % either normalize the weights (all positive)
        for n=nodes, gm.Alg.wt(n) = gm.Alg.wt(n) / wTot; end;
      else                                       % or set positive node to enforce unit sum
        gm.Alg.wt(nPos) = 1-wTot;
      end;
    end;  % (update weights)

  end;    % (if # nodes > 1)

  for j=1:nNodes,
    n = gm.Alg.nodes{i}(j);
    if (1) % TODO: (isempty(bel{j})),         % if we haven't pre-computed the belief,
      bel{j}=gm.Alg.theta{n};           %   compute it now (not including back msg)
      for c=gm.Alg.children{n}, bel{j}=bel{j}+gm.Alg.msgFwd{c}; end;
    else                                      % if we have, then just remove the back msg
      bel{j}=bel{j} - gm.Alg.msgBwd{n};
    end;                                      % and then eliminate X 
%if (n==481 || n==483 || n==1252 || n==1253 || n==859 || n==1115),
% fprintf('=====\n'); n, bel{j}, logsumexpPower(bel{j},X, 1.0/gm.Alg.nodes(n).wt), 
%end;
    gm.Alg.msgFwd{n} = logsumexpPower(bel{j},X, 1.0/gm.Alg.wt(n));
    bel{j}=log(factor());
    if (gm.Alg.parent(n) == 0),         % add roots to overall bound
      lnZ = lnZ + table(gm.Alg.msgFwd{n}); 
      %fprintf('Root %d => %f\n',n,table(gm.Alg.nodes(n).msgFwd));
    end; 
  end;

end; % (buckets / elimination order)
