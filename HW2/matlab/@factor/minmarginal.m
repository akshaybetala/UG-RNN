function f = minmarginal(F, vars)
% minmarginal: minimize the factor over all but a subset of its variables
% f = minmarginal(F,v) : compute the min-marginal of F wrt variables v, i.e., f = \min_{x\v} F(x)

% (c) Alexander Ihler 2010

f = min(F, vdiff(F.v,uint32(vars)));
