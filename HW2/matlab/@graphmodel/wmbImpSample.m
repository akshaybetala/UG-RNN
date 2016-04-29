function [ws, xs] = wmbImpSample(gm, nBatch)
% [w, x]=wmbImpSample(gm, nBatch) : get importance-weighted samples from a WMB
% Note: run wmbFwd(..) before using

if (~strcmp(gm.Alg.name,'WMB'))
  fprintf('Graphical model has not been specialized to the WMB algorithm\n');
  return;
end;

if (nargin < 2), nBatch = 1; end;

ws = zeros(nBatch, 1);        % save importance weights
xs = zeros(nBatch,gm.nVar);   % save sampled state vectors

for i=length(gm.Alg.order):-1:1,  % reverse elimination order
  X = gm.Alg.order(i); 
  nodes = gm.Alg.nodes{i};
  nNodes = length(gm.Alg.nodes{i});
  if (nNodes == 0) continue; end;
  qs = zeros(nBatch,nNodes);                         % store sampling probs
  nodewt = [gm.Alg.wt(nodes)];                 % get weights & compute
  cdf = [0;cumsum(nodewt(:))];                       %   weight distribution;
  [tmp, idx] = histc(rand(1,nBatch)*cdf(end),cdf);   % sample k from node idx(k)

  bel = cell(1,nNodes); 
  for j=1:nNodes,
    n = gm.Alg.nodes{i}(j);                   % compute belief at node n:
    bel{j}=gm.Alg.theta{n}; %+gm.Alg.nodes(n).msgBwd; % theta+bwd+\sum fwd
    for c=gm.Alg.children{n}, bel{j}=bel{j}+gm.Alg.msgFwd{c}; end;
    bel{j}=bel{j}-gm.Alg.msgFwd{n};                % remove forward msg (!)
    %bel{j}=bel{j}*(1.0/gm.Alg.nodes(n).wt);              % power 
    % should be (conditionally) normalized if msgFwd run prevously
    %exp(bel{j}*(1.0/gm.Alg.nodes(n).wt)),
  end;
  for k=1:nBatch,
    bk = bel{idx(k)};
    vk = vars(bk); vk=vk(vk~=X);
    pr1 = condition(bk, vk, xs(k,vk));
    pr = exp( pr1 * (1.0/gm.Alg.wt(nodes(idx(k)))) );
    xs(k,X) = sample(pr,1);
    for jj=1:nNodes,
      ws(k) = ws(k) + value(bel{jj},xs(k,vars(bel{jj})));   % numerator: f(x)
      % could condition but should be unneeded (TODO:check?)
      qs(k,jj) = value(bel{jj},xs(k,vars(bel{jj})))*(1.0/gm.Alg.wt(nodes(jj)));
    end;
  end;
  wv = [gm.Alg.wt(nodes)]';
  q = exp(qs) * wv;
  ws = ws - log(q+1e-300);   % denominator: q(x)
end; % (buckets / elimination order)
