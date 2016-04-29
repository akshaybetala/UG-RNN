function [X lnQ JG] = sample(gm, n, ord, ibound, JG)
% [X lnQ JG] = sample(gm, nSamples, order, ibound [, JG]) : draw samples from the graphmodel
% Arguments:
%   nSamples : # of samples to draw
%   order    : elimination order for gm (samples drawn in reverse order)
%   ibound   : ibound for approximate sampling, if gm does not admit exact inference
%   JG       : join graph data structure (overrides ibound if passed)
% Returns:
%   X        : samples
%   lnQ      : log probability of each sample drawn under the sampling distribution
%   JG       : join graph data structure
%
% TODO: Update to use new/better WMB structure?

if (nargin < 5) 
  [nil JG] = MBEstruct(gm,ord,ibound,'sum+');
end;

X   = zeros(n,max(JG.cliqElim));			% allocate storage for samples
lnQ = zeros(n,max(JG.cliqElim));			%   and their probability

for xv=ord(end:-1:1),   % in reverse elimination order, sample each variable
  R = min(find(JG.cliqElim == xv));		% use 1st region for sampling (?)

  F = factor();												% construct belief at that region
  for i=1:length(JG.factors{R}), F = F .* JG.factors{R}{i}; end;
  for e=find(JG.eDst==R), F = F .* JG.msg{e}; end;
  F = F ./ sum(F, xv);								% make a conditional distribution q(x|xPa)
  xPa = vars(F); xPa = xPa(xPa~=xv); 	% get parent variables, to condition on

  for i=1:n,													% sample a continuation of each tuple
    q = condition(F,xPa,X(i,xPa));		% 
    X(i,xv) = sample(q);
    lnQ(i,xv) = log(value(q,X(i,xv)));	% and save its probability
  end;
end;
lnQ = sum(lnQ,2);

