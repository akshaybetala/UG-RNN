"""
wmbe.py

Defines a weighted mini-bucket elimination procedure for approximate inference in MRFs

Version 0.0.1 (2015-09-28)
(c) 2015 Alexander Ihler under the FreeBSD license; see license.txt for details.
"""

import numpy as np
import time as time
from sortedcontainers import SortedSet;

from .factor import *
from .graphmodel import *


class wmbe(graphmodel, algorithm): 
    """A container class for a weighted mini-bucket based graphical model (jointree)"""

    def __init__(self,factors, elimOrder=None, elimOps=1):
        """Create the graphical model from a set of factors etc
        factors: list of factor objects comprising the model
        elimOrder: an elimination order or ordering method
        elimOp: elimination operator (0: max, or 1: sum) or list (one per variable)
        """
        self.order = []    # elimination order of the variables
        self.clique= []    # len-N list of (variable sets): clique of vars involved in node n
        self.theta = []    # len-N list of factors: current log-factor of node n
        self.belief= []    # len-N list of factors: current belief of node n
        self.parent= []    # len-N list of the parent of each node
        self.children = [] # len-N list of (list of children of this node)
        self.wt    = []  # length-N list of weight associated with each node
        self.msgFwd= []  # length-N list of forward message *from* this node to parent
        self.msgBwd= []  # length-N list of backward message *into* this node from parent
        #self.mb    = __mbNode(var, ...)?
        self.node  = []  # [ [nodes1] , [nodes2] , ... ] : B x Nb
        self.match = []  # [ [ [nodes1-match1], [nodes1-match2] ], [[nodes2-match1]], ... ] : B x M x Nm

        

"""
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
    gm.Alg.match{i}{end+1} = uint32([]);    % ensure at least one entry for consistency
  end;
  % (TODO) Remove non-unique match patterns...

  % (TODO) If desired, do moment matching on minibucket(i) list?

end;  % (elimination of Xi)
"""


def GibbsSampling( model, Query, state=None, nSamples=1 ):
  """'Functional' Gibbs sampling algo: """
  state = state if state is not None else [np.random.randint(Xi.states) for Xi in model.X]
  Query.table[:] = 0.0    # reset counts to zero
  for j in xrange(nSamples):
    for Xi in model.X:
      p = Factor([],1.0)
      for f in model.factorsWith(Xi):
        cvar = f.vars - [Xi]
        p *= f.condition2( cvar, [state[v] for v in cvar] )
      p /= p.sum()
      state[Xi] = p.sample()[0]
    Query[ tuple(state[v] for v in Query.vars) ] += 1.0
  Query /= nSamples
  return Query


class Query:
  """Defines a number of Monte Carlo "query" objects, for online estimation of various quantities.
     myQuery = Query.type(args)  :  create a new query of type "type"
     myQuery(x, w=1.0)           :  update query results after observing state "x" with optional weight w
     myQuery.value               :  access estimated query results

     Supported types:
       Query.marginals( GraphModel ) : estimate the marginals of all factors in a GraphModel object
       Query.expectation( G )        : estimate the expected value of G(x)
       Query.sequence( G )           : keep track of the entire trace of G(x) during the run
       Query.stateSequence()         : keep track of the trace of x during the run
  """
  # TODO: Query.queryList?  Query.maximum( F )?

  @staticmethod
  def marginals(model):
    def track(x,w=1.0):
      # update all factors: *(step-1)/step then [x] + 1.0/step
      for f in track.factors:
        f *= float(track.step-1)/track.step
        f[tuple(x[v] for v in f.vars)] += 1.0/track.step
      track.step += 1
    # initialize object and return it
    track.factors = [Factor(f.vars) for f in model.factors] # init factors
    track.step = 1
    return track

  @staticmethod
  def expectation(G):
    def track(x,w=1.0):
      track.value *= float(track.step-1)/track.step
      track.value += w*G(x)/track.step
      track.step += 1
    # initialize object and return it
    track.value = 0.0 
    track.step = 1
    return track

  @staticmethod
  def sequence(G):
    def track(x,w=1.0):
      track.w.append(w)
      track.values.append(G(x))
    # initialize object and return it
    track.w = []
    track.values = []
    return track

  @staticmethod
  def stateSequence(): return Query.sequence(lambda x: x)

#  @staticmethod
#  def stateSequence():
#    def track(x,w=1.0):
#      track.w.append(w)
#      track.x.append(x)
#    # initialize object and return it
#    track.w = []
#    track.x = []
#    return track


  

#  use functions?  "marginals(model)(x,w) => update list of marginals in model"
#      "mapconfig(F)(x,w) => eval F(x) & update if larger"
#      "expectation(G)(x,w) => update expectation of G(x)   (needs step #?)
#      "stateSequence()(x,w) => append x to state sequence list
#      "expectSequence(G)(x,w) => append G(x) to sequence

class Gibbs:
  """Gibbs sampling algorithm"""

  state = []
  EQ = []

  def __init__(self, model=None, state=None):
    """Create a gibbs sampling approximate inference object"""
    self.model = model   # save reference to graphical model
    self.state = state if state is not None else [np.random.randint(Xi.states) for Xi in model.X]
    #if (state is None and model is not None): x = [np.random.randint(Xi.states) for Xi in model.X]

# thoughts...
#  query only single variable marginals?
#  query all cliques in the model?
#  query one or some cliques in the model?
#  query expectation of some function?
#  query "trace" of a function (possibly identity => state sequence)
#  query generic operator, e.g., keep argmax_xi f(xi) ?
#     (could use for expectations with callable object?)



  def setQueryCliques(self, cliqueList ):
    self.queryCliques = cliqueList
    # reset results?

  def setQueryFunctions(self, functionList ):
    self.queryFunctions = functionList
    # reset results?

  def setTraceFunctions(self, functionList ):
    self.traceFunctions = functionList

  def reset(self):
    step = 0
    for fa in self.queryCliques: fa *= 0.0
    

  def run(self, maxIter=-1, maxTime=float('inf')):
    if (maxIter == -1) and (maxTime==float('inf')): raise ValueError('No maximum iterations or time set')
    import itertools
    start_time = time.time()
    for j in itertools.count(1):
      if j > maxIter: break;
      if time.time() - start_time > maxTime: break;
      # take a pass through each variable and sample its value given the others:
      for Xi in self.model.X:
        p = Factor([],1.0)
        for f in self.model.factorsWith(Xi):
          cvar = f.vars - [Xi]
          p *= f.condition( cvar, [self.state[v] for v in cvar] )
        p /= p.sum()
        self.state[Xi] = p.sample()[0]
      # after each pass, evalute desired queries:
      # print self.state
      for fa in self.queryCliques:
        fa *= float(j-1)/j
        fa[ tuple(self.state[v] for v in fa.vars) ] += 1.0/j



class Metropolis:
  """Metropolis-Hastings algorithm"""

  # init: need model (or just model eval f'n?), proposal draw & eval, queries

  # classfunction default proposal? randomly choose k variables & flip them?

  # 





### Forward / rejection sampling  (for BNs only?)
###    Verify GM is a BN; sample x; reject if not Xe; update query

### Likelihood weighting (BNs only?  Or any with topo order conditionals?)
###    Verify GM is a BN; sample Xh & eval Xe

### Importance sampling (generic form?)
###   f(x), q(x) draw and eval; queries
###   





 
