function f = marginal(F, vars)
% marginal: sum out all but a subset of variables in the factor
% f = marginal(F,v) : compute the marginal of F wrt variables v, i.e., f = \sum_{x\v} F(x)

% (c) Alexander Ihler 2010

f = sum(F, vdiff(F.v,uint32(vars)));
