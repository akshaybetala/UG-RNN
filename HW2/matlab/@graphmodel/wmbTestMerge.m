function [mergeScores, mergeCliques] = wmbTestMerge(gm, iBound,sBound)

if (nargin<3) sBound=inf; end;
if (nargin<2) iBound=inf; end;

mergeCliques = {};
mergeScores = [];
logDims = log2(gm.dims);

for i=1:length(gm.Alg.order),
X = gm.Alg.order(i);
nodes = gm.Alg.minibucket(i).nodes;

  for j=1:length(nodes),
    for k=j+1:length(nodes), 
      n1 = nodes(j); n2=nodes(k);
      cliqueNew = vunion( gm.Alg.nodes(n1).clique, gm.Alg.nodes(n2).clique );
      if (length(cliqueNew) > iBound+1), continue; end;     % if ibound violated, skip
      if (sum(logDims(cliqueNew)) > sBound), continue; end; % if sbound violated, skip
      if (length(cliqueNew) == max(length(gm.Alg.nodes(n1).clique),length(gm.Alg.nodes(n2).clique)))
        fprintf('yes?\n'); pause;
        mergeScores(end+1) = inf;
        mergeCliques{end+1} = cliqueNew;
        continue;
      end;
      msgNew = gm.Alg.nodes(n1).theta + gm.Alg.nodes(n2).theta;
      for c=gm.Alg.nodes(n1).children, msgNew=msgNew+gm.Alg.nodes(c).msgFwd; end;
      for c=gm.Alg.nodes(n2).children, msgNew=msgNew+gm.Alg.nodes(c).msgFwd; end;
      wtNew = gm.Alg.nodes(n1).wt + gm.Alg.nodes(n2).wt;

      %% max+ bound: optimum of combined function vs separate optima of individuals:
      score = max(msgNew + gm.Alg.nodes(n1).msgBwd - gm.Alg.nodes(n2).msgFwd) - max(gm.Alg.nodes(n1).msgFwd + gm.Alg.nodes(n1).msgBwd);
      score = score * 2.^(length(logDims)-i);
      score = (length(logDims)-i);

if (0)
      msgNew = logsumexpPower(msgNew,X, 1.0/wtNew);
      %% max+ bound: optimum of combined function vs separate optima of individuals:


      %score = max(msgNew)-max(gm.Alg.nodes(n1).msgFwd)-max(gm.Alg.nodes(n2).msgFwd);
      % sum+ bound (optimistic): lnZ of combined function vs summed lnZ of individuals:
      score = logsumexp(msgNew)-logsumexp(gm.Alg.nodes(n1).msgFwd)-logsumexp(gm.Alg.nodes(n2).msgFwd);
      % sum+ bound (pessimistic): sum^0 = max of updated clique, with reparameterization:
      score = max(msgNew-gm.Alg.nodes(n1).msgFwd-gm.Alg.nodes(n2).msgFwd);

%      fprintf('Considering merge %d+%d: ', n1,n2);
%      fprintf('%d\n',score);
end;

      mergeScores(end+1) = score;
      mergeCliques{end+1} = cliqueNew;
    end;
  end;
end;

