"""
montecarlo.py

Defines several Monte Carlo and MCMC routines for approximate inference in graphical models

Version 0.0.1 (2015-09-28)
(c) 2015 Alexander Ihler under the FreeBSD license; see license.txt for details.
"""

import numpy as np
import time as time
from sortedcontainers import SortedSet;

from .factor import *
from .graphmodel import *


####### Basic sampling #########################
#
# sample( factors ) -- in sequence, draw all vars not yet sampled, conditioned on sampled values; also return q(x)
# 
# importance sample:  sample(fs) => x,q ; value(x) => p  => x, w=p/q
#
# rejection sample: sample(fx) => x,q ; reject if r > p/q*R  where q*R > p for all x   (shortcut?)
#
# gibbs sample: 
#
#
# fwdSample (?) - in a BN, sample each factor in sequence; weight by evidence  => = IS w/ conditionals
#

#
# Basic sampling function interface:
#   x,logp = sample();    # return configuration "x" and log p(x), the log-probability of drawing that sample
# For Markov chain Monte Carlo:
#   xb,logp = sample(xa); # return sampled transition xa -> xb and log p(xa->xb), the probability of that transition 
#

# Basic query interface:
#   Q = Query( f )                   # to compute expectation of f(x): E_p[f] 
#   Q = Query( [f1, f2, ...] )       # to compute E_p[fi] for each fi in the list
#   Q = QueryMarginals( factorlist ) # to compute p(x_a) for each factor f_a(x_a) in the list
#
# A query is an object Q=Query(..) with:
#   Q() or Q[i]   : return the current estimate, or ith estimate from a list
#   Q.update(x,w) : update the expectations' estimates by observing state x with weight w
#   Q.wvar()      : return the empirical variance of the weights encountered (optional)


class Query:
    """Defines a Monte Carlo "query" object, for estimation of various quantities from sequences of (weighted) states.
       Q = Query( f )   # function f(x), estimate expectation E_p[ f(x) ]
       Q = Query( [f1,f2,...] )   # functions fi(x), estimate each expectation

       An object with a query interface should include at least:
           Q()  : return the current estimate(s)
           Q[i] : return the ith estimate (if a list of estimates)
           Q.update(x,w) : update the estimates after observing state "x" with weight "w"
           Q.reset()     : reset / re-initialize estimates
       It may also have
           Q.nsamples  : total number of samples (calls to update)
           Q.wtot      : total of weights seen during calls to update
           Q.neff      : number of "effective" samples
           Q.wvar      : variance of weight values
    """
    def __init__(self, functions):
        self.isList = True
        if (not hasattr(functions,"__getitem")): 
            functions = [functions]
            self.islist = False
        self.functions = functions
        self.sums  = [0.0] * len(functions)
        self.nsamples = 0.0
        self.wtot  = 0.0  # save weight total
        self.w2tot = 0.0  # save weight^2 total
        # TODO: should probably be in log-weight domain

    def update(self,x,w):
        for i,f in enumerate(self.functions):
            self.sums[i] += f(x)
        self.nsamples += 1
        self.wtot += w
        self.w2tot+= w**2

    def __getitem__(self,i):
        return self.sums[i]/self.wtot

    def __call__(self):
        to_return = [ s / self.wtot for s in self.sums ]
        if not self.isList: to_return = to_return[0]
        return to_return

    @property
    def wvar(self):
        return (self.w2tot - self.wtot**2)/self.nsamples

    @property
    def neff(self):
        return self.wtot**2 / self.w2tot

###################################################### 
class QueryMarginals:
    """Specialized Monte Carlo "query" object for marginal probabilities of factors
       Q = QueryMarginals( factorlist ) # estimate marginal p(x_a) for each factor f_a(x_a) in list
    """
    def __init__(self, factorlist):
        self.marginals = [ Factor(f.vars,0.0) for f in factorlist ]
        self.nsamples = 0.0
        self.wtot = 0.0
        self.w2tot= 0.0

    def update(self,x,w):
        for mu in self.marginals:
            mu[ tuple(x[v] for v in mu.vars) ] += w
        self.nsamples += 1
        self.wtot += w
        self.w2tot += w**2

    def __getitem__(self,i):
        return self.marginals[i]/self.wtot

    def __call__(self):
        return [ mu / self.wtot for mu in self.marginals ]

    # TODO: should probably just inherit from Query
    @property
    def wvar(self):
        return (self.w2tot - self.wtot**2)/self.nsamples

    @property
    def neff(self):
        return self.wtot**2 / self.w2tot




def GibbsSampling( model, query, state=None, stopSamples=1, stopTime = inf ):
    """Gibbs sampling procedure for graphical model "model" with query object "query"
    """
    state = state if state is not None else [np.random.randint(Xi.states) for Xi in model.X]
    # TODO: timing check
    for j in xrange(nSamples):
        # TODO if out of time, break
        for Xi in model.X:
            p = Factor([],1.0)
            for f in model.factorsWith(Xi):
                cvar = f.vars - [Xi]
                p *= f.condition2( cvar, [state[v] for v in cvar] )
            p /= p.sum()
            state[Xi] = p.sample()[0]
        query.update(state, 1.0)
    return query






"""
def GibbsSampling( model, Query, state=None, nSamples=1 ):
  ""'Functional' Gibbs sampling algo: 
  ""
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
  ""Defines a number of Monte Carlo "query" objects, for online estimation of various quantities.
     myQuery = Query.type(args)  :  create a new query of type "type"
     myQuery(x, w=1.0)           :  update query results after observing state "x" with optional weight w
     myQuery.value               :  access estimated query results

     Supported types:
       Query.marginals( GraphModel ) : estimate the marginals of all factors in a GraphModel object
       Query.expectation( G )        : estimate the expected value of G(x)
       Query.sequence( G )           : keep track of the entire trace of G(x) during the run
       Query.stateSequence()         : keep track of the trace of x during the run
  ""
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

"""
  
#  use functions?  "marginals(model)(x,w) => update list of marginals in model"
#      "mapconfig(F)(x,w) => eval F(x) & update if larger"
#      "expectation(G)(x,w) => update expectation of G(x)   (needs step #?)
#      "stateSequence()(x,w) => append x to state sequence list
#      "expectSequence(G)(x,w) => append G(x) to sequence

class Gibbs:
  ""Gibbs sampling algorithm""

  state = []
  EQ = []

  def __init__(self, model=None, state=None):
    ""Create a gibbs sampling approximate inference object""
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
          p *= f.condition2( cvar, [self.state[v] for v in cvar] )
        p /= p.sum()
        self.state[Xi] = p.sample()[0]
      # after each pass, evalute desired queries:
      # print self.state
      for fa in self.queryCliques:
        fa *= float(j-1)/j
        fa[ tuple(self.state[v] for v in fa.vars) ] += 1.0/j



class Metropolis:
  ""Metropolis-Hastings algorithm""

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





 
