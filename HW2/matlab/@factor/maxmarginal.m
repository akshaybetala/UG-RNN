function f = maxmarginal(F, vars)
% maxmarginal: maximize over all but a subset of variables of a factor
% f = maxmarginal(F,v) : compute the max-marginal of F wrt variables v, i.e., f = \max_{x\v} F(x)

% (c) Alexander Ihler 2010

f = max(F, vdiff(F.v,uint32(vars)));
