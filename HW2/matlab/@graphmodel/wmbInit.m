function [gm lnZ] = wmbInit(gm, order,iBound,sBound,elimOp, options)
% [gm lnZ] = wmbInit(gm, order, iBound, sBound, elimOp [,options])
%  Create & initalize data structure for weighted mini-bucket elimination
%  TODO: more header comments

if (~isempty(gm.Alg.name) && ~strcmp(gm.Alg.name,'WMB')), 
  error('Graphical model already specialized to a different algorithm?') 
end;

gm.Alg.name = 'WMB';                   % specialize graphical model to WMB algorithm:
gm.Alg.order = uint32(order);          %   save the elimination order (must be fixed)
gm.Alg.clique{1}   = uint32([]);       %   list of nodes: clique variables,
gm.Alg.theta{1}    = log(factor()); %  (log) potential function 
gm.Alg.belief{1}   = log(factor()); %  (log) potential function 
gm.Alg.parent(1)   = uint32(0);  %     parent node (send forward message to)
gm.Alg.children{1} = uint32([]); %     children nodes (rec'v fwd messages from)
gm.Alg.wt(1)       = 1;          %     elimination weight
gm.Alg.msgFwd{1} = log(factor());%     forward messages from node n to its parent
gm.Alg.msgBwd{1} = log(factor());%     backward msg *into* node n from its parent
gm.Alg.minibucket(1).var   = order(1); %   minibucket structure: variable eliminated,
gm.Alg.nodes{1} = uint32([]);    %     nodes corresponding to MB partitioning
gm.Alg.match{1} = {{uint32([])}};%     node sets on which to perform matchings

REFILL = false; % TODO: option parsing

n = 1;                                 % n: current/next node position to add
nVar = gm.nVar;                        % # of variables in the model
logDim = log2(gm.dims);                % log # of states per variable
lnZ = 0.0;                             % initialize bound

factors = gm.factors;                  % first convert factors to log-factors 
for f=1:length(factors), factors{f}=log(factors{f}); end;
gmo = graphmodel(factors);             % Copy the original factors for manipulation, and
fIsMsg = zeros(1,length(gmo.factors)); %   keep track of whether each factor is a message, 
fSrc   = 1:length(gmo.factors);        % and its "source" (node/orig factor #)
% TODO: worry about vacant positions?

for i=1:length(order),
  X = order(i);
  gm.Alg.minibucket(i).var = X;
  gm.Alg.nodes{i}=uint32([]);
  gm.Alg.match{i}={};

  % Get bucket for current variable:
  bucketIds = withVariable(gmo, X);    % get current factors that contain X
  bucket = getFactor(gmo, bucketIds);  %   (their positions, the factors, and # of factors)  
  nFact = length(bucket); if (nFact==1) bucket={bucket}; end;
  if (nFact==0) continue; end;         % if no factors involve X, just skip it (!!!)
  gmo=removeFactor(gmo,bucketIds);     % remove those factors from the collection
  bIsMsg = fIsMsg(bucketIds);          %   and keep track of their sources & "ismsg"
  bSrc   = fSrc(bucketIds);

  % Do moment matching among factors before partitioning? (for value-based partitions?)
  % (???)

  % Select allocations into buckets
  %   keep track of "sources" of partition: orig factor # & msg origination
  % TODO: Easy way? 
  nBucket = length(bucket);
  sz=zeros(1,nBucket); for j=1:nBucket, sz(j)=nvar(bucket{j}); end;  % (???) or use table size?
  [nil srt] = sort(sz,'descend');                      % sort factors in terms of 
  bucket=bucket(srt); bIsMsg=bIsMsg(srt); bSrc=bSrc(srt);  % decreasing size
  for c=1:nBucket,
    clique{c}=uint32([]);  % start with empty partitions
    msgIn{c} =uint32([]);  % no messages in yet
    factIn{c}=uint32([]);  % nor factors assigned yet
    factTo{c}=uint32([]);  % (reverse clique/factor identifiers)
  end;
  for j=1:nBucket,         % for each factor,
    for c=1:nBucket,       %   find a clique to put it in:
      vnew = vunion( vars(bucket{j}), clique{c} );
      logSz= sum(logDim(vnew)); 
      iBoundTmp = max([iBound,length(clique{c})-1,length(vars(bucket{j}))-1]);    % can't be smaller than
      sBoundTmp=max([sBound,sum(logDim(clique{c})),sum(logDim(vars(bucket{j})))]);% current clique or factor
      if ((length(vnew) <= iBoundTmp + 1) && logSz <= sBoundTmp ),
        factTo{j} = uint32(c);             % assign factor j to clique c
        clique{c} = vnew;                  %  increasing the clique size
        if (bIsMsg(j)) msgIn{c} =[msgIn{c}  bSrc(j)];  % which msgs & 
        else           factIn{c}=[factIn{c} bSrc(j)];  %  factors are in c?
        end;
        break;             % (exit from finding clique for j)
      end;
    end;
  end;
  nClique=nBucket; for c=1:nBucket, if (isempty(clique{c})), nClique=c-1; break; end; end;
  % TODO: check factTo{} to make sure they all ended up somewhere...

  % Refill: run through factors again in reverse order & ask if can be added
  if (REFILL)
    for j=nBucket:-1:1,
      for c=1:nClique,
        vnew = vunion( vars(bucket{j}), clique{c} );
        logSz= sum(logDim(vnew));
        iBoundTmp = max(iBound,length(clique{c})-1); sBoundTmp=max(sBound,sum(logDim(clique{c})));
        if ((length(vnew) <= iBoundTmp + 1) && logSz <= sBoundTmp ),  % if we fit,
          factTo{j} = vunion(factTo{j}, uint32(c)); %  keep track of that,
          clique{c} = new;                         %  and increase the clique size
        end;
      end;
    end;
  end;  % (if refill)

  % Assign: allocate factors into clique potentials
  theta = cell(1,nClique);  bel = cell(1,nClique); % create empty potential functions 
  for c=1:nClique, theta{c} = log(factor()); end;  %   for the model parameters and beliefs
  bel = theta;
  for j=1:nBucket,                                 % "simple assignment": place factor &
    c = factTo{j}(1);                              %   message into first assigned clique
    bel{c} = bel{c}+bucket{j};
    if (~bIsMsg(j)) theta{c}=theta{c}+bucket{j}; end;
  end;
  % TODO: FIX
  %for j=1:nBucket,                                % "uniform assignment": place factor & msg
  %  ftmp = bucket{j} ./ length(factTo{j});         % assign each factor to their cliques equally
  %  for c=factTo{j}, theta{c}=theta{c}+ftmp; end;  % in equal proportion
  %end;

  % For each new node n, update its data structure and add result msg to gmo:
  for c=1:nClique,
    % n = next node position in list (constructed in order)
    gm.Alg.wt(n) = 1.0/nClique;            % (TODO) for sum+ (upper bound) only (!!!)
    if (strcmp(elimOp,'max+')),                  % (TODO): max+ version:
      gm.Alg.wt(n) = 1e-6; % eps;          % (TODO): eps? numerical stability, etc?
    end;
    if (strcmp(elimOp,'sum-') && nClique>1),      % (TODO): sum- version:
      if (c==1) gm.Alg.wt(n) = 2.0; else gm.Alg.wt(n)=-1.0/(nClique-1); end; 
    end;
    gm.Alg.clique{n} = clique{c};          % store clique variables
    gm.Alg.theta{n}  = theta{c};           %   clique parameters, & forward msg
    %if (any(clique{c}~=vars(bel{c}))), clique{c}, vars(bel{c}), pause; end;
    gm.Alg.msgFwd{n} = logsumexpPower(bel{c},X, 1.0/gm.Alg.wt(n)); % Eliminate
    gm.Alg.parent(n) = 0;                  % no parent (yet)
    gm.Alg.msgBwd{n} = log(factor());      % no backward msg (yet)
    % (TODO)  For assigned factors, note f->n map & node's list
    gm.Alg.children{n} = msgIn{c};         % get children (nodes sending fwd message)
    for cc=msgIn{c}, gm.Alg.parent(cc)=n; end; % and fix their parent pointer
    gm.Alg.nodes{i} = [gm.Alg.nodes{i} n]; % add node to MB struct

    if (isempty(vars(gm.Alg.msgFwd{n}))),  % if it's a scalar outcome,
      lnZ = lnZ+table(gm.Alg.msgFwd{n});   %  add it to the bound (=> no parent)
    else                                         % else, add new factor to incremental graph
      [gmo, pos] = addFactor(gmo, gm.Alg.msgFwd{n});
      fIsMsg(pos)=1; fSrc(pos)=n;                  % keeping track of the source info
    end;
    n = n+1;                                     % and advance to next node position
  end; % (creation of each node)

  % (TODO) Generate match structure from factors & "factTo"
  for j=1:nBucket, % any factor,
    if (length(factTo{j})>1)
      gm.Alg.match{i}{end+1} = gm.Alg.nodes{i}(factTo{j});
    end;
  end; 
  if (nClique>1)
    gm.Alg.match{i}{end+1} = uint32(gm.Alg.nodes{i}); % and "all"
	else 
		gm.Alg.match{i}{end+1} = uint32([]);		% ensure at least one entry for consistency
  end;
  % (TODO) Remove non-unique match patterns...

  % (TODO) If desired, do moment matching on minibucket(i) list?

end;  % (elimination of Xi)


