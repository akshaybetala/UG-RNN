function [f,alpha]=normalize(f,method)
% construct a normalized factor under one of several interpretations
% [f,alpha] = normalize(F, method) : normalize F according to some metric.  method can be:
%  'sum' : ensure f sums to one  (e.g., for probabilities): F=alpha*f;
%  'max' : make the largest value of f equal to zero  (e.g., for logarithmic values) F=f+alpha

% (c) Alexander Ihler 2010

  if (nargin<2) method='sum'; end;
  switch (lower(method))
    case 'sum', 
      alpha=sum(f.t(:)); 
      if (alpha<=0) error('factor:normalize: factor is non-positive'); end;
      f.t = f.t./alpha;
    case 'max', alpha=max(f.t(:)); f.t = f.t-alpha;
  otherwise, error('factor:normalize: unknown normalization method');
  end;

