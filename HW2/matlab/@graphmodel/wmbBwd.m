function gm = wmbBwd(gm)
% gm=wmbBwd(gm [,options]) : compute backward messages for weighted mini-bucket elimination

if (~strcmp(gm.Alg.name,'WMB'))
  fprintf('Graphical model has not been specialized to the WMB algorithm\n');
  return;
end;


for i=length(gm.Alg.minibucket):-1:1,  % for each bucket (= collection of minibuckets)
  X = gm.Alg.order(i); %minibucket(i).var;     % variable to be eliminated
  nodes = gm.Alg.nodes{i};
  nNodes = length(gm.Alg.nodes{i});  % # of minibucket partitions in bucket
  bel = cell(1,nNodes);             % allocate storage for beliefs

  if (length(gm.Alg.nodes{i})>1) % if more than one mini-bucket partition:

    if (1)                          % (TODO) if matching, compute beliefs jointly
      for j=1:nNodes,
        n = gm.Alg.nodes{i}(j);                   % compute belief at node n:
        bel{j}=gm.Alg.theta{n}+gm.Alg.msgBwd{n}; % theta+bwd+\sum fwd
        for c=gm.Alg.children{n}, bel{j}=bel{j}+gm.Alg.msgFwd{c}; end;
        bel{j}=bel{j}*(1.0/gm.Alg.wt(n));              % and power (TODO?)
        bel{j}=bel{j} - logsumexp(bel{j});                   % and normalize
      end;

      % (TODO) Can do reparameterization & weight updates here....

    end;  % (if 1)
  end;    % (if # nodes > 1)

  % Compute backward messages from each node n to its children c
  for j=1:nNodes,
    n = gm.Alg.nodes{i}(j);
    wt_n = gm.Alg.wt(n);
    if (1) %isempty(bel{j})),            % (TODO) if no belief computed, compute now 
      bel{j}=gm.Alg.theta{n}+gm.Alg.msgBwd{n}; % theta+bwd+\sum fwd
      for c=gm.Alg.children{n}, bel{j}=bel{j}+gm.Alg.msgFwd{c}; end;
      bel{j}=bel{j} - logsumexp(bel{j});                   % and normalize
    end;
    for c=gm.Alg.children{n},
      wt_c = gm.Alg.wt(c);
      vElim = vdiff( gm.Alg.clique{n} , gm.Alg.clique{c} );

      %msg = logsumexp( bel{j}*(1.0/wt_n) + (gm.Alg.nodes(c).msgFwd*(-1.0/wt_c)) ,vElim)*wt_c;
      % (TODO) : numerical stability issues??? treatment of exact zeros?
      %mu = min( wt_c , wt_n );  
      %msg = bel{j}*(mu/wt_n) + (gm.Alg.nodes(c).msgFwd*(-mu/wt_c));
      %msg = msg - max(msg);
      %msg = logsumexp( msg*(1.0/mu) ,vElim )*wt_c;
      % ALT:
      msg = bel{j} - max(bel{j}); % TODO: just do this earlier
      msg = logsumexp( msg*(1.0/wt_n) ,vElim )*wt_c - gm.Alg.msgFwd{c};
 
      gm.Alg.msgBwd{c} = msg;
    end;                            
    bel{j}=log(factor());           % clear out belief
  end;  % (backward msgs)


end; % (buckets / elimination order)




