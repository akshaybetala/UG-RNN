function [pi,cliques,ptparent] = order(gm,method,priority)
% [pi,cliq,parent] = order(model,method,priority) : heuristic selection of variable elimination order
%  method = 'MinFill','WtMinFill'    (min-fill and wieghted min-fill heuristic)
%           'MinWidth','WtMinWidth'  (min-width and wieghted min-width heuristic)
%           'Random'                 (random elimination ordering)
%  priority = 1xnVar list of "priority" of elimination -- all nodes with priority 1 1st, then 2, etc.
%             (useful for restricted orderings for mixed elimination problems; default = all 1)
%  outputs:
%    pi     : the ordering, from first eliminated to last eliminated
%    cliq   : the cliques induced by the elimination process
%    parent : pseudo-tree parent; the clique that this clique is eliminated into (zero if a root)
%

%nVar = max(gm.variables);
nVar = gm.nVar;
dims = gm.dims; dims(dims==0)=1;
nFactors = length(gm.factors);
fixNbrs = false;
scores=zeros(1,nVar);                        % create vector of scores for each varible
if (nargout >= 2) cliques=cell(1,nVar); end;  % create container for cliques if requested
if (nargin < 3) priority = ones(1,nVar); end; % all variables have equal priority unless specified

switch(lower(method)),
  case 'minfill',    score = @(a,b,c) scoreMinFill(a,b,c,0,dims);  fixNbrs=true;
  case 'wtminfill',  score = @(a,b,c) scoreMinFill(a,b,c,1,dims);  fixNbrs=true;
  case 'minwidth',   score = @(a,b,c) scoreMinWidth(a,b,c,0,dims); fixNbrs=false;
  case 'wtminwidth', score = @(a,b,c) scoreMinWidth(a,b,c,1,dims); fixNbrs=false;
  case 'random',     score = @(a,b,c) a; scores=rand(size(scores)); fixNbrs=false;
  otherwise, error('Unknown ordering method');
end;

% Construct adjacency lists for each variable
adj = cell(1,nVar); for i=1:nVar, adj{i}=uint32(i); end;
pi = zeros(1,nVar,'uint32');    % placeholder for the order
for j=1:nFactors, if (~isempty(gm.factors{j})),  % go through each factor and figure
  v=variables(gm.factors{j});                    %   out the adjacency of the MRF
  for i=v,
    adj{i} = vunion(adj{i},v);
  end;
end; end;

% Compute initial heuristic values
scores=score(scores,adj,1:nVar);

% Now pick the order:
for ii=1:nVar,                               %  
  minPri = min(priority);
  scoresAllow = scores; scoresAllow(priority > minPri) = inf;
  [tmp,p]=min(scoresAllow);                  % choose node with min heuristic value
  %p = find(scores==tmp); p=p(ceil(length(p)*rand));  % select at random from among equal-value nodes
  priority(p)=inf;
  p=uint32(p); pi(ii)=p;  
  v = adj{p}; fix=uint32([]);                % update connectivity of neighbors
  if (nargout > 1) cliques{ii}=v; end;       % and list of cliques constructed
  for i=v,
    adj{i} = vdiff(vunion(adj{i},v),p);      % fix up adjacency after eliminating
    if (fixNbrs) fix=vunion(fix,adj{i}); end;  % and add neighbor & nbrs' nbrs to the fix-up list
  end;
  if (~fixNbrs) fix = v; end;
  scores = score(scores,adj,fix);
  scores(p) = inf;           % mark the chosen node & don't reselect
end;


if (nargout >= 3)   % requested pseudotree information
  rpi(pi)=1:nVar;
  ptparent = zeros(1,nVar,'uint32');
  for ii=1:nVar,
    tmp = rpi(vdiff(cliques{ii},pi(ii)));            % find ordering of neighbors at pi(ii)'s elimination
    if (~isempty(tmp)) ptparent(ii)=min(tmp); end;   % leave zero if no parent
  end;
end;



function scores = scoreMinFill(scores,adj,fix,weighted,dims)
  scores(fix) = 0;
  for i=fix,
    for ii=adj{i},
      if (weighted) scores(i)=scores(i) + prod(dims(vdiff(adj{i},adj{ii}))); 
      else          scores(i)=scores(i) + length(vdiff(adj{i},adj{ii}))/2;
      end;
    end; 
  end;

function scores = scoreMinWidth(scores,adj,fix,weighted,dims)
  for i=fix,
    if (weighted) scores(i)=prod(dims(adj{i})); 
    else          scores(i)=length(adj{i});
    end;
  end;



