function [bel, del] = wmbBelief(gm, X)
% [gm, lnBound]=wmbFwd(gm [,options]) : compute forward messages for weighted mini-bucket elimination

if (~strcmp(gm.Alg.name,'WMB'))
  fprintf('Graphical model has not been specialized to the WMB algorithm\n');
  return;
end;

X = uint32(X);
i = find(gm.Alg.order == X);
del = factor([],0.0);
if (isempty(i) || length(gm.Alg.nodes{i})==0), bel=factor([],1.0); return; end;

  n = gm.Alg.nodes{i}(1);                   % compute belief at node n:
  bel=gm.Alg.theta{n}+gm.Alg.msgBwd{n}; % theta+bwd+\sum fwd
  for c=gm.Alg.children{n}, bel=bel+gm.Alg.msgFwd{c}; end;
  bel=bel*(1.0/gm.Alg.wt(n));              % and power (TODO?)
  bel=logsumexp(bel, vdiff(vars(bel), X));
  bel=bel - logsumexp(bel);             % normalize
  bel = exp(bel);

for j = 1:length(gm.Alg.nodes{i}),
  n = gm.Alg.nodes{i}(j);                   % compute belief at node n:
  delta{j}=gm.Alg.theta{n}-gm.Alg.msgFwd{n}; % theta+bwd+\sum fwd
  for c=gm.Alg.children{n}, delta{j}=delta{j}+gm.Alg.msgFwd{c}; end;
  delta{j}=max(delta{j}, vdiff(vars(delta{j}), X));
  del = del + delta{j};
end


