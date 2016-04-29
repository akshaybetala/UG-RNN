function gm = trainCD(gm, data, nIter, T, initStep)
% gm = trainCD(gm, data, nIter, T, initStep) : train graphical model with contrastive divergence
%
% TODO: comment, check, fix bugs, etc (everything!)

if (nargin < 5) initStep = .1; end;
if (nargin < 4) T = 1; end;
if (nargin < 3) nIter = 100; end;

F0 = cell(size(gm.factors));
for f=1:length(gm.factors)
  v = vars(gm.factors{f});
  F0{f} = empirical(gm.factors{f},data(:,v));
end;

samp=data;
for iter=1:nIter,
  step = initStep/iter;
  for i=1:size(data,1),
    samp(i,:) = gibbs(gm, 1, T-1, 1, data(i,:));
    fprintf('.');
  end;
  fprintf('\n');
  for f=1:length(gm.factors),
    v = vars(gm.factors{f});
    Ft{f} = empirical(gm.factors{f},samp(:,v));
    gm.factors{f} = gm.factors{f} * exp(step*(F0{f}-Ft{f}));
    gm.factors{f} = normalize(gm.factors{f});
  end;
end;



