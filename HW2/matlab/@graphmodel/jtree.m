function [lnZ, bel] = jtree(gm,ord)
% [lnZ, bel] = jtree(gm [,order]) -- compute partition function and marginals exactly using junction tree
%  bel : marginal probabilities for each factor in the graphical model
%  lnZ : log partition function of the graphical model

nFactors = length(gm.factors);
gmCopy = gm;                                % (save a copy)
[ord,cliques,parent]=order(gm,'MinFill');   % get an elimination order, associated cliques and jt parents
cliquefns = cell(size(cliques));            % storage for clique beliefs
msgs = cell(size(cliques));                 % storage for downward messages
children = cell(size(cliques));             % junction tree children (computed from parent)
map = zeros(1,nFactors);                    % where is factor i in the junction tree?
for i=1:length(cliquefns), cliquefns{i}=factor([],1.0); children{i} = uint32([]); end;
lnZ = 0.0;

for i=1:length(ord),                        % run through variable elimination ordering
  X = ord(i);                               % 
  fids = withVariable(gm,X);                % getting factors that go with this clique
  map(fids) = i;                            %   and saving a map from those factors to this clique
  for j=1:length(fids), cliquefns{i} = cliquefns{i}*getFactor(gm,fids(j)); end;
  gm = removeFactor(gm,fids);               %
end;

for i=1:length(ord),                        % message passing down the tree
  p=parent(i);
  Z = sum(cliquefns{i});                    % normalize belief collected so far
  lnZ = lnZ + log(Z);                       %   and save normalization constant
  cliquefns{i} = cliquefns{i}/Z;   
  if (p~=0)
    msgs{i} = marginal(cliquefns{i},cliques{p});  % compute downward message
    cliquefns{p} = cliquefns{p} * msgs{i};        % and place into parent clique f'n
    children{p} = vunion(children{p},uint32(i));  % (& save child ids)
  end;
end;

for p=length(ord):-1:1,                     % message passing back up the tree
  Z = sum(cliquefns{p});                    %
  lnZ = lnZ + log(Z);                       % re-normalize belief & save norm constant
  cliquefns{p} = cliquefns{p}/Z;
  for i=children{p},                        % compute upward messages (exclude downward msg)
    msg = marginal(cliquefns{p}/msgs{i},cliques{i});
    cliquefns{i} = cliquefns{i} * msg;      % and place into child clique f'ns
  end;
end;
  
bel=cell(1,nFactors);                       % compute marginals of each original factor in model
for i=1:nFactors, if (~isempty(getFactor(gmCopy,i))),
  bel{i} = marginal(cliquefns{map(i)}, vars(getFactor(gmCopy,i)));
end; end;


