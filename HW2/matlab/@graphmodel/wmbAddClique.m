function gm = wmbAddClique(gm, clique)
% gm=wmbAddClique(gm, clique) : add a clique to the WMB approximation
% 

if (~strcmp(gm.Alg.name,'WMB'))
  fprintf('Graphical model has not been specialized to the WMB algorithm\n');
  return;
end;

[nil, priority] = sort(gm.Alg.order);    % find out which variable will be eliminated first
i = min(priority(clique)); 
[porder, corder] = sort(priority(clique));%  and the order in which they are eliminated

n = length(gm.Alg.nodes)+1;  % (TODO) use vacancy stack...
nDone = [];

% walk through the cliques that will be created by eliminating the added one:
for l=1:length(clique),
  X = clique(corder(l));                 % which bucket will each one go in?
  i = porder(l); % = priority(X)         %  "" ""
  subsumed = 0;                          %   once there's already a mini-bucket that
  for j=gm.Alg.minibucket(i).nodes,      %   contains it, we can stop adding cliques
    if (isempty(vdiff(clique(corder(l:end)),gm.Alg.nodes(j).clique))),
      subsumed=j; break;
    end; 
  end;
  if (subsumed~=0),                      % if we found one, point last created node into it
    if (~isempty(nDone)),
      gm.Alg.nodes(j).children(end+1) = nDone(end); 
      gm.Alg.nodes(nDone(end)).parent = j;
    end;
    break; 
  end;
  % Otherwise, we need to create a new clique:
  gm.Alg.nodes(n).clique = sort(clique(corder(l:end)));
  gm.Alg.nodes(n).theta  = log(factor());
  gm.Alg.nodes(n).belief = log(factor());
  gm.Alg.nodes(n).parent = 0;               % no parent (for now)
  gm.Alg.nodes(n).children = [];          
  if (~isempty(nDone)), 
    gm.Alg.nodes(n).children = nDone(end);  % child is just previously created node (for now)
    gm.Alg.nodes(nDone(end)).parent = n;    %  and this is its parent 
  end;
  gm.Alg.nodes(n).wt = 1e-6; %eps;          % (TODO): initialize zero weight (increases if merge later)
  gm.Alg.nodes(n).msgFwd = log(factor());   % initialize: blank fwd/bwd msgs
  gm.Alg.nodes(n).msgBwd = log(factor());
  nDone(end+1) = n;                         % append this to the list of created nodes
  n = length(gm.Alg.nodes)+1;               % (TODO) use vacancy stack
end;

for l=1:length(nDone),
  X = clique(corder(l));                 % which bucket will each one go in?
  i = porder(l); % = priority(X)         %  "" ""
  n = nDone(l);
  for j=gm.Alg.minibucket(i).nodes(end:-1:1), % check if any nodes are subsumed by the new one
    if (isempty(vdiff(gm.Alg.nodes(j).clique,gm.Alg.nodes(n).clique))),
      % if j subsumed (=> remove j), j's children become n's children
      gm.Alg.nodes(n).children = [gm.Alg.nodes(n).children gm.Alg.nodes(j).children];
      for c=gm.Alg.nodes(j).children, gm.Alg.nodes(c).parent = n; end;
      p=gm.Alg.nodes(j).parent;               % remove j from parent's child list
      if (p~=0), 
        gm.Alg.nodes(p).children = gm.Alg.nodes(p).children(gm.Alg.nodes(p).children~=j);
      end;
      %fprintf('Removing %d with parent %d, children [',j,p); 
      %for c=gm.Alg.nodes(j).children, fprintf('%d,',c); end; fprintf(']\n');
      % absorb j's parameters, weight into node n
      gm.Alg.nodes(n).theta = gm.Alg.nodes(n).theta + gm.Alg.nodes(j).theta;
      gm.Alg.nodes(n).wt    = gm.Alg.nodes(n).wt + gm.Alg.nodes(j).wt;  
      % initialize forward message using j's (will be replaced in fwd pass)
      gm.Alg.nodes(n).msgFwd = gm.Alg.nodes(n).msgFwd + gm.Alg.nodes(j).msgFwd;
      % (TODO) initialize msgBwd?
      % If j's parent will not be subsumed, need to reparameterize to maintain current sol'n:
      %   n's parent -= msgFwd ; j's parent += msgFwd 
      %   (TODO): should also be sure to "match" between { pa(n) & pa(j) }
      if (p~=0),
        pn=gm.Alg.nodes(n).parent;
        gm.Alg.nodes(pn).theta = gm.Alg.nodes(pn).theta - gm.Alg.nodes(j).msgFwd;
        gm.Alg.nodes(p).theta  = gm.Alg.nodes(p).theta  + gm.Alg.nodes(j).msgFwd;
      end;
      % remove j from the minibucket list (nodes in which Xi is eliminated)
      gm.Alg.minibucket(i).nodes(gm.Alg.minibucket(i).nodes==j)=[]; 
      % remove j from list of matchings (replacing with n) & ensure that n doesn't appear twice
      for m=1:length(gm.Alg.minibucket(i).match),
        gm.Alg.minibucket(i).match{m}( gm.Alg.minibucket(i).match{m}==j ) = n;  % j -> n
        gm.Alg.minibucket(i).match{m} = unique(gm.Alg.minibucket(i).match{m});  % make unique
      end;
    end;
  end;
  gm.Alg.minibucket(i).nodes(end+1) = n;   % add n to the minibucket list for Xi
end;

% (TODO): Need to update match structure 
%  (1) remove singleton and duplicate match patterns
%  (2) order in some canonical way?
%  (3) add new match patterns: (a) reparameterization, (b) any new messages?, (c) all?
% (TODO): Reparameterization change for child-not-contained case (+ ensure match added)

