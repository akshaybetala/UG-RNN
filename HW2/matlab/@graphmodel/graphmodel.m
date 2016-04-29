function gm = graphmodel(flist);
%    graphmodel({factors})  : construct from a cell array of factors
%
  if (nargin == 0) flist={}; end;

  % create empty graph model -- enforces correct order
  %                     % Representing the distribution itself:
  gm.factors={};        %  - list of factors (probability tables)
  gm.vacant=uint32([]); %  - empty entries of factors{:}  (FILO queue)
  gm.nVacant=0;         %  - how many vacant entries do we have? (queue size)
  gm.nVar=0;            %  - # of variables (actually largest var index)
  gm.dims=[];           %  - dimensions of each variable (sized by Nv? or variables?)
  gm.vNbrs=[];          %  - factors that depend on each variable
  gm.Alg.name='';       % Algorithmic-specific storage
  %                     % Representing relationships among components
  gm.eVacant = uint32(4:-1:1);
  gm.neVacant = 4;
  gm.edges = zeros(4,4,'uint32');
  gm.adj = {};
  %fg.nbrs=[];            %  - factors that neighbor on one another (edges of graph)
  %fg.edgeIndex=[];      %  - lookup table for edge-based entries
  %fg.eIdx=sparse([],[],[]);  % (i->j) => e
  %fg.eSrc=uint32([]);        % e => i
  %fg.eDst=uint32([]);        % e => j
  %
  %fg.msg=[];            % For message passing algorithms
  %fg.msgNew=[];          % 
  %fg.beliefs=[];        % 
  %fg.wts=[];            % 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Initialize factor graph : factors and variable neighborhoods
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Nf = length(flist);

  gm.factors=flist; gm.vacant=zeros(1,16,'uint32'); %zeros(size(flist),'uint32');
  nVar=0;
  for i=1:length(flist), 
    if (~isempty(flist{i})),
      % need to ensure copy operator on factor list? !!!
      var=variables(flist{i});
      %gm.variables=vunion(gm.variables,var);
      nVar=max([nVar,max(var)]);
    else
      gm.factors{i}=factor();   % different (zero) if log-factors... (?)
      gm.nVacant=gm.nVacant+1;
      gm.vacant(gm.nVacant)=i;
    end;
  end;
  gm.nVar=nVar;
  %nVar = length(gm.variables);

  vNbrs=cell(1,nVar); dim=zeros(1,nVar); %fEmpty=uint32([]);
  for i=1:length(flist)
    if (~isempty(flist{i})),
      var=variables(flist{i});
      dim(var)=dims(flist{i});
      %vidx = find(vmember(gm.variables,var));
      %dim(vidx)=dims(flist{i});
      %for j=vidx, vNbrs{j}=vunion(vNbrs{j},uint32(i)); end;
      for j=var, vNbrs{j}=vunion(vNbrs{j},uint32(i)); end;
    end;
  end;
  gm.vNbrs=vNbrs; gm.dims=dim; %gm.fEmpty=fEmpty;

  gm.adj = cell(1,length(gm.factors));
  for i=1:length(gm.factors), gm.adj{i}=uint32([]); end;

  gm=class(gm,'graphmodel');

%  gm=consolidate(gm);
